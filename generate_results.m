clc;close all; clear;
% General Simulation Parameters
params.Ps =    100e-3;  % 100 milliwatt
params.simulation_area_side = [-500 500]; % square of side 1 km
params.space_realizations = 1000;
params.time_slots = 1000  ;

% Network generation
disp('Multiple association in Ultra-dense network, ICC 2016');
disp('generating network nodes ...');
simulation_area = (params.simulation_area_side(2) - params.simulation_area_side(1)).^2; 
disp(['simulation area ', num2str(simulation_area/1e6) , ' km^2']);
disp('generating the channel: Rayleigh fading');
params.H = exprnd(1,100.5e6,1); % Rayleigh fading channels

% Probability of Idle state results
%params.k =[0.1 0.5 1 2:2:8 10:5:35 40:10:90  100:20:200] ;
%params.la_u =  0.00003;  % 300 users/km^2
params.N = [1 5 10];  % number of RF chains of the UE = number of small cells in a virtual cell

%  disp('generating idle mode probability results...');
%  results.Po = generate_idle_mode_probability_results(params);

% AER results
params.la_s =[0.001 0.005  0.01:0.01:0.09 0.1]; 
params.la_u =  0.0003;  % 300 users/km^2
params.simulation_area_side = [-500 500]; % square of side 1 km
params.space_realizations = 100;
params.time_slots = 100  ;
params.N = 3 ;
params.alpha = 4;
params.rho =   0.5e-6;
disp('generating average ergodic rate results...');
results.AER = generate_AER_results_mod(params);


% params.la_s = 0.01;
% params.rho_dBm = -60:5:-30;
% params.rho = 10.^(params.rho_dBm/10) * 1e-3;
% params.simulation_area_side = [-1000 1000]; % square of side 1 km
% params.space_realizations = 1000;
% params.time_slots = 100  ;
% disp('generating average ergodic rate results (rho)...');
% results.AER_rho = generate_AER_rho_results(params);

% ASE results
%params.la_s =[0.001 0.005  0.01:0.03:0.09 0.1:0.04:0.3 0.4 0.5]; 
% params.la_s =[0.001 0.005 0.01 0.05 0.1 0.5]; 
% params.la_u =  0.0003;  % 300 users/km^2
%  params.simulation_area_side = [-500 500]; % square of side 1 km
%  params.space_realizations = 10;
%  params.time_slots = 10  ;
% % 
% % %params.N = 1;
% % params.N = [1 2 5];
%  params.alpha = 4;
%  params.rho =   1e-7;
% % %params.backhaul_capacity = 100e6;
% % %params.bandwidth = 20e6;
% % results.ASE = generate_ASE_results(params);
% 
% params.N = [1 2 5];
% params.la_s = 0.01;
% params.la_u =  0.0001:0.0003:0.001; 
% results.ASE_lau = generate_ASE__lau_results(params);
% %Coverage results
% % params.la_s =[0.001:0.002:0.009 0.01:0.02:0.09 0.1:0.2:0.9 1] 
% % params.la_u =  0.0003;  % 300 users/km^2
% % params.simulation_area_side = [-100 100]; % square of side 1 km
% % params.space_realizations = 10;
% % params.time_slots = 10  ;
% % disp('generating coverage results...');
% % generate_coverage_results(params);
