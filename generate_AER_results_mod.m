function result = generate_AER_results_mod(params)
la_s = params.la_s ;
exact = zeros(numel(params.N),numel(la_s));
simul = zeros(numel(params.N),numel(la_s));
M= params.N;
for j = 1:M
    for i = 1:numel(la_s)
        fprintf('\n')
        disp(['Nth nearest: ' , num2str(j) , '  la_s: ', num2str(la_s(i)) ]);
        exact(j,i) = comp_rate(j,la_s(i),params);
        simul(j,i) = sim_rate(j,la_s(i),params);
    end
end

simul = simul / log(2);
exact = exact / log(2);

result = {exact , simul} ;
h = plot(la_s,exact,'-r',la_s,simul,'ko');
%h = plot(la_s,exact,'-r',la_s,simul,'ko');
set(h,'MarkerSize',10);
set(h,'LineWidth',4);
xlabel('\textbf{$\lambda_s$}','interpreter','latex');
ylabel('Average spectral efficiency ($bps/Hz$)','interpreter','latex');
title('Average spectral efficiency ($\alpha = 4 , \sigma^2 = 0$ )','interpreter','latex');
legend(h([1 6]),{'Exact','Simulation'},'FontSize',14,'FontWeight','bold');
set(gca, 'FontSize', 20);
set(gca, 'FontWeight', 'Bold');

end

function rate = comp_rate(j,la_s,params)
alpha = params.alpha;  %4;
la_u = params.la_u ; % 0.0003;
rho = params.rho; %1e-6;  % 0.1 micro watt
Ps = params.Ps;%100e-3;  % 100 milliwatt
% idle mode probability calculation
%la_s = 0.1;

k = la_s/la_u;
%M = params.N;
p = (3.5*k / (params.N + 3.5*k))^3.5;
%p = 1- 3/k;
% neighborhood radius calculations
a = (rho/Ps)^(-1/alpha);
%N_avg = pi*la_s * a.^2;
% noise term
%.* exp(-(sig^2 / Ps) * r.^alpha .* (exp(t) - 1) )
%j = 1;
f_j =  @(t,r)  2/gamma(j)*(pi*la_s)^j * r.^(2*j-1) .* exp(- pi*la_s * r.^2 .* (1 + (1 - p) * sqrt(exp(t) - 1) .* (pi/2 - atan((a^2 ./ r.^2) .* 1./sqrt(exp(t) - 1)))))   ;
rate = integral2(f_j,0,inf,0,inf);
end


function rate = sim_rate(j,la_s,params)

%alpha = params.alpha;  %4;
la_u = params.la_u ; % 0.0003;
rho = params.rho; %1e-6;  % 0.1 micro watt
Ps = params.Ps;%100e-3;  % 100 milliwatt
%runs = params.space_realizations * params.time_slots;
simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
%simulation_area =  pi * params.simulation_radius^2;
mu = la_s * simulation_area;
nu = la_u * simulation_area;
rates = zeros(params.space_realizations,params.time_slots);
count = 0;
m = 1;
H = params.H; % Rayleigh fading channels


while (m <= params.space_realizations)
    if(mod(m,params.space_realizations/100) == 0)
        fprintf('|');
    end
    N_cells = poissrnd(mu);
    N_users = poissrnd(nu);
    R = zeros(N_users,N_cells);  % received signal matrix
    A = zeros(N_users,N_cells);  % association matrix
    S = zeros(N_users,1);        % received signal from server for each user
    I = zeros(N_users,1);        % interference signal from interferers for each user
    %Server = zeros(N_users,params.N);
    Ncu = zeros(N_users,1);
    cell_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
    user_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users, 2);
    % Compute distances to the typical user
    %r = cell_pos(:,1).^2 + cell_pos(:,2).^2;
    %##### Edge Effect correction
    %         a = params.simulation_area_side(2) - params.simulation_area_side(1);
    %         r = min(cell_pos(:,1),a - cell_pos(:,1)).^2 + min(cell_pos(:,2),a - cell_pos(:,2)).^2;
    r = pdist2(cell_pos,user_pos,'euclidean').^2;
    for t = 1:params.time_slots
        h = reshape(H(t:N_cells*N_users+t-1),N_cells,N_users);
        R = Ps * h .* (1./(r.*r));   % special optimization in case of alpha = 4
        neighbours  = getStrongest(R,j);
        interferers = R(R<rho);
        % SIR computations
        %k = la_s/la_u;
        %p = (3.5*k / (1 + 3.5*k))^3.5; % probability of idle mode
        %p = (3.5*k / (M + 3.5*k))^3.5;
        try
            S = neighbours(j);
        catch
            m;
            t;
            count = count + 1;
            continue;
            
        end
        I = (1 - p) * sum(interferers);
        %I = (1 - p)^2 * sum(interferers);
        SIR = S/I;
        
        % Rate computations
        rates(m,t) = log(1 + SIR) ;
    end
    m = m + 1;
end

rate = sum(sum(rates)) / numel(rates);
end

function [signals , servers ] = getStrongest(neighbors, N)
N=2;
N_user = size(neighbors,2);
signals = zeros(N_user,N);
servers = zeros(N_user,N);
for i = 1:N
    ne = neighbors;
    [signals(:,i) , ind] = max(neighbors);
    neighbors(ind,i) = 0;
    servers(:,i) = ind;
end
end

