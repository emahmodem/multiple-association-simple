function result = generate_AER_rho_results(params)
    rho = params.rho ;
    M = params.N;
    exact = zeros(M,numel(rho));
    simul = zeros(M,numel(rho));

    for j = 1:M
        for i = 1:numel(rho)
            fprintf('\n')
            disp(['Nth nearest: ' , num2str(j) , '  rho: ', num2str(rho(i)) ]);
            exact(j,i) = comp_rate_rho(j,rho(i),params);
            simul(j,i) = sim_rate_rho(j,rho(i),params);
        end
    end

    simul = simul / log(2);
    exact = exact / log(2);

    result = {exact , simul} ;
    h = plot(params.rho_dBm,exact,'-r',params.rho_dBm,simul,'ko');
    set(h,'MarkerSize',10);
    set(h,'LineWidth',4);
    xlabel('\textbf{$\rho$}','interpreter','latex');
    ylabel('Average spectral efficiency ($bps/Hz$)','interpreter','latex');
    title('Average spectral efficiency ($\alpha = 4 , \sigma^2 = 0$ )','interpreter','latex');
    legend(h([1 6]),{'Exact','Simulation'},'FontSize',14,'FontWeight','bold');
    set(gca, 'FontSize', 30);
    set(gca, 'FontWeight', 'Bold');

end
 

function servers = getStrongest(neighbors, N)
    servers = zeros(N,1);
    for i = 1:N
        [servers(i) , ind] = max(neighbors);
        neighbors(ind) = 0;
    end
end

function rate = comp_rate_rho(j,rho,params)
M = params.N;
alpha = params.alpha;  %4;
la_u = params.la_u ; % 0.0003;
la_s = params.la_s; %1e-6;  % 0.1 micro watt
Ps = params.Ps;%100e-3;  % 100 milliwatt
% idle mode probability calculation 
%la_s = 0.1;

k = la_s/la_u;
p = (3.5*k / (M + 3.5*k))^3.5;

% neighborhood radius calculations
a = (rho/Ps)^(-1/alpha);
%N_avg = pi*la_s * a.^2;
% noise term
%.* exp(-(sig^2 / Ps) * r.^alpha .* (exp(t) - 1) )
%j = 1;
f_j =  @(t,r)  2/gamma(j)*(pi*la_s)^j * r.^(2*j-1) .* exp(- pi*la_s * r.^2 .* (1 + (1 - p) * sqrt(exp(t) - 1) .* (pi/2 - atan((a^2 ./ r.^2) .* 1./sqrt(exp(t) - 1)))))   ;
rate = integral2(f_j,0,inf,0,inf);

end


function rate = sim_rate_rho(j,rho,params)
    M = params.N;
    %alpha = params.alpha;  %4;
    la_u = params.la_u ; % 0.0003;
    la_s = params.la_s;
    %rho = params.rho; %1e-6;  % 0.1 micro watt
    Ps = params.Ps;%100e-3;  % 100 milliwatt
    %runs = params.space_realizations * params.time_slots;
    simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
    %simulation_area =  pi * params.simulation_radius^2;
    mu = la_s * simulation_area;
    rates = zeros(params.space_realizations,params.time_slots); 
    count = 0;
    m = 1;
    H = params.H; % Rayleigh fading channels
    while (m <= params.space_realizations)
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        N_cells = poissrnd(mu);
        cell_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        % Compute distances to the typical user
        r = cell_pos(:,1).^2 + cell_pos(:,2).^2;

        for t = 1:params.time_slots
            h = H(t:N_cells+t-1);
            R = Ps * h .* (1./(r.*r));   % special optimization in case of alpha = 4
            %neighbours = sort(R(R>=rho), 'descend');
            neighbours = getStrongest(R,params.N);
            interferers = R(R<rho);
            % SIR computations
            k = la_s/la_u;
            p = (3.5*k / (M + 3.5*k))^3.5; % probability of idle mode
            try
                S = neighbours(j);
            catch  
            m
            t
            count = count + 1
            continue;

            end
            I = (1 - p) * sum(interferers);
            %I = (1 - p)^2 * sum(interferers);
            SIR = S/I;
        
            % Rate computations
            rates(m,t) = log(1 + SIR) ; 
        end
        m = m + 1;
        %progressbar(m/params.space_realizations);
    end
    
      rate = sum(sum(rates)) / numel(rates);
      % rate = sum(rates) / numel(rates);
end


