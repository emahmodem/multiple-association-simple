function result = generate_idle_mode_probability_results(params)
    simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
    nu = params.la_u * simulation_area;
    la_u = params.la_u;
    la_s = params.k * la_u ; 
    
    simul = zeros(numel(params.N),numel(la_s));
    for n = 1: numel(params.N)
        for l = 1:numel(la_s)
            disp(['N: ' , num2str(params.N(n)) , '  la_s: ', num2str(la_s(l)) , ' ' , num2str(round((n*l)/(numel(la_s) * numel(params.N) ) * 100)) , '%']);
            disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
            disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
            probs = zeros(params.space_realizations,1);
            mu = la_s(l) * simulation_area;
            for m = 1:params.space_realizations
                if(mod(m,params.space_realizations/100) == 0)
                    fprintf('|');
                end
                N_cells = poissrnd(mu);
                N_users = poissrnd(nu);
                A = zeros(N_cells,1); % association vestor; 1: cell is active , 0: cell is idle
                cell_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
                user_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users, 2);
                distances = pdist2(user_pos,cell_pos);
                for u = 1:N_users
                    for a = 1:params.N(n);
                        [~,i] = min(distances(u,:));
                        A(i) = 1;
                        distances(u,i) = inf;  
                    end
                end
                probs(m) = 1 - mean(A);
            end
            simul(n,l) = mean(probs); 
            fprintf('\n')
        end
    end
    k = la_s ./ la_u;
    exact = zeros(numel(params.N),numel(k));
    approx = zeros(numel(params.N),numel(k));
    for n = 1:numel(params.N)
        exact(n,:) = ((3.5 * k) ./ (params.N(n) + 3.5 * k)).^ 3.5;
        approx(n,:) = 1 - params.N(n)./k;
    end
    result = {exact , simul , approx} ;
    h = plot(k,exact,'-r',k,simul,'ro',k,approx,'k--');
    set(h,'MarkerSize',15);
    set(h,'LineWidth',4);
    xlabel('\textbf{$\kappa = \frac{\lambda_s}{\lambda_u}$}','interpreter','latex');
    ylabel('p_o');
    title('');
    legend(h([1 4]),{'Exact','Simulation'},'FontSize',14,'FontWeight','bold');
    set(gca, 'FontSize', 20);
    set(gca, 'FontWeight', 'Bold');
    ylim([0 1]);
end
     