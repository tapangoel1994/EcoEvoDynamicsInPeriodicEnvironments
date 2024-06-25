%% This script generates the parameters (life history traits and simulation parameters) that are kept constant across simulations.

%% Life history parameters (units of hours, micrograms and mL). 
params.J = 0; %ug/mL-h
params.conversion_efficiency = 5e-7; %ug/cell
params.d_R = 0; % per hour
params.mu_max = 1.2; % per hour
params.R_in = 4; %ug/mL
params.alpha_l = 0;
params.alpha_e = 0;
params.alpha_i = 0;

params.d_S = .2; %per hour
params.d_E = .2; %per hour
params.d_L = .2; %per hour
params.d_I = .2; %per hour
params.m = 1/24; %per hour

params.phi = 3.4e-10; %mL/hr
params.lambda = 2; %per hour
params.eta = 1; %per hour
params.bet = 50;


%% Simulation parameters:
params.flask_volume = 500; %volume in mL
params.dt = 1/30; % hours
MaxCycles = 100;
steadystatethresh = 1e-1/params.flask_volume; % concentration difference below which two concentrations are treated as identical

%% Numerical method related parameters
options = odeset('AbsTol',1e-8,'RelTol',1e-8,'NonNegative',1:10); %Options for the ODE function call