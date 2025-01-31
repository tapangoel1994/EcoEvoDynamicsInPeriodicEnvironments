%% This script generates the parameters (life history traits and simulation parameters) that are kept constant across simulations.

%% Life history parameters (units of hours, micrograms and mL). 
params.J = 0; %ug/mL-h %resource inflow rate
params.conversion_efficiency = 5e-7; %ug/cell %% amount of resources that need to be consumed to create 1 cell
params.d_R = 0; % per hour %decay rate of resources
params.mu_max = 1.2; % per hour %maximal cell growth rate
params.R_in = 4; %ug/mL %resource density at half max cell growth rate
params.alpha_l = 0; % fitness advantage of lysogens over susceptible cells
params.alpha_e = 0; % fitness advantage of exposed cells over susceptible cells
params.alpha_i = 0; % fitness advantage of lytically infected cells over susceptible cells

params.d_S = .04; %per hour %death rate of susceptible cells
params.d_E = .04; %per hour %death rate of exposed cells
params.d_L = .04; %per hour %death rate of lysogens
params.d_I = .04; %per hour %death rate of lytically infected cells
params.m = 1/24; %per hour %virion decay rate

params.phi = 3.4e-10; %mL/hr %adsorption rate
params.lambda = 2; %per hour %transition rate (from exposed to lytic or lysogenic cell)
params.eta = 1; %per hour %lysis rate
params.R_eta = params.R_in; %reducing in lysis rate
params.bet = 50; % burst size


%% Simulation parameters:

params.dt = 1/30; % hours
MaxCycles = 50000; % Max number of cycles to steady state before while loop terminates
InvasionCycles = 10;% Number of cycles in each set to evaluate transients during invasion
criticaldensitythreshold = 1e-3; % concentration difference below which two concentrations are treated as identical in per mL
params.flask_volume = 1/criticaldensitythreshold; %volume of system in mL

%% Numerical method related parameters
options = odeset('AbsTol',1e-8,'RelTol',1e-8,'NonNegative',1:10); %Options for the ODE function call

%% Initial conditions
R0 = 1e2; %ug/mL
S0 = 1e7; % cells/mL
Va_0 = 1e4; %virions/mL
Vb_0 = 10*criticaldensitythreshold; %viral genome copies/mL

%% Transfer parameters
q_R = 0; % fraction of resources transferred from one cycle to next
q_S = 0; % fraction of susceptible cells transferred from one cycle to next
q_E = 0; % fraction of exposed cells transferred from one cycle to next
q_I = 0; % fraction of lytically infected cells transferred from one cycle to next

PGamma = struct('lytic',[0 0],'temperate',[.92 .0063],'lysogenic',[1 0]); %temperate strategy is optimal for T = 24hr, q_L = 0.2, q_V = 0
Strategies = {'lytic','temperate','lysogenic'};
%save('fixedparameters.mat');