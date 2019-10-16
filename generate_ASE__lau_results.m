function result = generate_ASE__lau_results(params)
    la_u = params.la_u ;
    exact = zeros(numel(params.N),numel(la_u));
    simul = zeros(numel(params.N),numel(la_u));
    for m = 1:numel(params.N)
        M = params.N(m); 
        for i = 1:numel(la_u)
            fprintf('\n')
            disp(['M = : ' , num2str(M) , '  la_u: ', num2str(la_u(i)) ]);
            exact(m,i) = comp_ASE(M,la_u(i),params);
            simul(m,i) = sim_ASE(M,la_u(i),params);
        end
        
    end
    simul = simul / log(2);
    exact = exact / log(2);

    result = {exact , simul} ;
    h = plot(la_u*1e6,exact,'-r',la_u*1e6,simul,'ko');
    grid on;
    set(h,'MarkerSize',10);
    set(h,'LineWidth',4);
    xlabel('\textbf{$\lambda_u$}','interpreter','latex');
    ylabel('Area spectral efficiency ($bps/Hz/km^2$)','interpreter','latex');
    title('Area spectral efficiency ($\alpha = 4 , \sigma^2 = 0$ )','interpreter','latex');
    legend(h([1 6]),{'Exact','Simulation'},'FontSize',20,'FontWeight','bold');
    set(gca, 'FontSize', 30);
    set(gca, 'FontWeight', 'Bold');

end
 
function ASE = comp_ASE(M,la_u,params)
    alpha = params.alpha;  %4;
    la_s = params.la_s ; % 0.0003;
    rho = params.rho; %1e-6;  % 0.1 micro watt
    Ps = params.Ps;%100e-3;  % 100 milliwatt
    % idle mode probability calculation 
    %la_s = 0.1;

    k = la_s/la_u;
    p = (3.5*k / (M + 3.5*k))^3.5;
    %p = 1- 3/k;
    % neighborhood radius calculations
    a = (rho/Ps)^(-1/alpha);
    %N_avg = pi*la_s * a.^2;
    % noise term
    %.* exp(-(sig^2 / Ps) * r.^alpha .* (exp(t) - 1) )
    %j = 1;
    rate = zeros(1,M);
    for j = 1:M
        f_j =  @(t,r)  2/gamma(j)*(pi*la_s)^j * r.^(2*j-1) .* exp(- pi*la_s * r.^2 .* (1 + (1 - p) * sqrt(exp(t) - 1) .* (pi/2 - atan((a^2 ./ r.^2) .* 1./sqrt(exp(t) - 1)))))   ;
        %f_j =  @(t,r)  2/gamma(j)*(pi*la_s)^j * r.^(2*j-1) .* exp(- pi*la_s * r.^2 .* (1 + (1 - p) * sqrt(exp(t) - 1) .* (pi/2 - atan(                 1./sqrt(exp(t) - 1)))))   ;
        rate(j) = integral2(f_j,0,inf,0,inf);
    end
    
    % backhaul limitation
    %rate = min(rate .* params.bandwidth , params.backhaul_capacity) ./ params.bandwidth;
    ASE = (la_u * 1e6) * sum(rate);
    
end


function ASE = sim_ASE(M,la_u,params)
    %alpha = params.alpha;  %4;
    la_s = params.la_s ; % 0.0003;
    rho = params.rho; %1e-6;  % 0.1 micro watt
    Ps = params.Ps;%100e-3;  % 100 milliwatt
    %runs = params.space_realizations * params.time_slots;
    simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
    %simulation_area =  pi * params.simulation_radius^2;
    mu = la_s * simulation_area;
    nu = la_u * simulation_area;
    rates = zeros(params.space_realizations,params.time_slots);
    
    H = params.H; % Rayleigh fading channels

    for m = 1:params.space_realizations
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        N_cells = poissrnd(mu);
        cell_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        N_users = poissrnd(nu);
        user_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users, 2);
        association = zeros(N_users,N_cells);
        user_rates = zeros(N_users,M);
            for t = 1:params.time_slots
                for u = 1:N_users
                    r = (user_pos(u,1) - cell_pos(:,1)).^2 + (user_pos(u,2) - cell_pos(:,2)).^2;
                    h = H(t:N_cells+t-1);
                    R = Ps * h .* (1./(r.*r));   % special optimization in case of alpha = 4
                    [neighbours , ind] = getStrongest(R,M);
                    association(u,ind) = 1;
                    %rho = neighbours(M);
                    interferers = R(R<rho);
                    % SIR computations
                    k = la_s/la_u;
                    %p = (3.5*k / (1 + 3.5*k))^3.5; % probability of idle mode
                    p = (3.5*k / (M + 3.5*k))^3.5;
                    S = neighbours;
                    I = (1 - p) * sum(interferers);
                    %I = (1 - p)^2 * sum(interferers);
                    SIR = S./I;

                    % Rate computations
                    for a = 1:M
                        user_rates(u,a) = log(1 + SIR(a)) / sum(association(u,:));
                    end
                end
                rates(m,t) = sum(sum(user_rates)); 
             end 
                
                
    end
         ASE = (sum(sum(rates)) / numel(rates)) / (simulation_area *1e-6);       
end
        
    
    
      


function [servers,ind] = getStrongest(neighbors, N)
    servers = zeros(N,1);
    ind = zeros(N,1);
    for i = 1:N
        [servers(i) , ind(i)] = max(neighbors);
        neighbors(ind(i)) = 0;
    end
end

