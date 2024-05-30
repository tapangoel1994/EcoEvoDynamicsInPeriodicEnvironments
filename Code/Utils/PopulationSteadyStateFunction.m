%%Code takes one-host one-virus system and sweeps over (q,gamma) keeping
%%other parameters fixed. Uses an RSEILV model with no MOI dependence. 
% Obtain appropriate steady state densities for
%%every (q,gamma) pair and generate a heatmap. 

%% Date Created: 1/23/2024
%% Author: Tapan Goel

%%Inputs:
% CyclePeriod - Duration of Single Cycle in hours.
% p_L - fraction of lysogens being passaged between cycles.
% p_V - fraction of virions being passaged between cycles.
% Gamma - array of Gamma values at which steady states are to be evaluated.
% Q - array of Q values at which steady states are to be evaluated. 
% numNodes - number of nodes being used for parallel computing. Do not
% exceed 32 (depending on cluster being used).

%%Output:
% Script generates a .mat file named:
% "CyclePeriod=<CyclePeriod>,S0=<InitialHostDensity>,V0=<InitialViralDensity>,p_L=<p_L>,p_V=<p_V>.mat"
% that contains all the variables in the workspace from the simulation.
% Note that the output file has variables for the steady state values of
% each cell type for all strategies evaluated but does not have the
% dynamics for the strategies.

function PopulationSteadyStateFunction(CyclePeriod,p_L,p_V,Gamma,Q,numNodes)

addpath('Utils');

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
params.q = [0 0];
params.gamma = [0 0];

%% simulation parameters:
params.flask_volume = 500; %volume in mL
params.dt = 1/30; % hours
params.T = CyclePeriod; % hours
params.t_vals = transpose(0:params.dt:params.T); % time
MaxCycles = 50000;

%% filter parameters
p_S = 0;
p_E = 0;
p_I = 0;
p_L = p_L;
p_V = p_V;
TransferMatrix = diag([0 p_S p_E p_E p_I p_I p_L p_L p_V p_V]);

%% Numerical method related parameters
options = odeset('AbsTol',1e-8,'RelTol',1e-8,'NonNegative',1:10); %Options for the ODE function call
steadystatethresh = 1e-1/params.flask_volume; % concentration difference below which two concentrations are treated as identical

SteadyStateDensityTemp = zeros(length(Q)*length(Gamma),10);
SteadyStateDensity = zeros(length(Q),length(Gamma),10);
SSCycles = zeros(length(Q)*length(Gamma),3); %% Variable to store q,gamma,#cyclestoSteadyState


%% initial conditions
R0 = 1e2; %initial resource amount in ug/mL ( 500 mL flask)
S0 = 1e7; %Initial concentration of susceptibles in flask (per mL)
V01= 1e4; %initial concentration of virus in flask (per mL)
V02 = 0;

%% Initiate parallel pool
poolobj = parpool(numNodes);

tic
parfor ii = 1:length(Q)*length(Gamma)
        Params = params;
        [j,i]=ind2sub([length(Gamma),length(Q)],ii);
        Params.q = [Q(i) 0];
        Params.gamma = [Gamma(j) 0];
        
        %% First cycle
        x0 = [R0 S0 zeros(1,6) V01 V02];
        [t_vals, y] = ode113(@ODE_RSEILV_2Species, Params.t_vals, x0, options, Params);
        
        %TimeSeries = y;
        %PlotTimeSeries_2species(t_vals,y);
        %drawnow;
        
        %% Cycles till steady state or till MaxCycles
        
        iter = 1;
        steadyrep = 0;
        
        x0 = [R0 S0 zeros(1,8)] + y(end,:)*TransferMatrix;
        
        %% loop continues till:
        %1. The initial population for the new epoch becomes (almost) equal to the
        %initial population for the previous epoch for atleast 10 cycles OR,
        %2. The total number of cycles hits the max number of cycles.
        
        while ( (sum(abs(x0 - y(1,:)) > steadystatethresh) >= 1) || steadyrep < 10) && iter < MaxCycles+10
            
            if((sum(abs(x0 - y(1,:)) > steadystatethresh) < 1))
                        steadyrep = steadyrep+1;
            end
        
        iter = iter+1;
        
        [t_vals, y] = ode113(@ODE_RSEILV_2Species, Params.t_vals, x0, options, Params);
        x0 = [R0 S0 zeros(1,8)] + y(end,:)*TransferMatrix;
        %TimeSeries = [TimeSeries;y];
        %PlotTimeSeries_2species(t_vals,y);
        %drawnow;
        end
        SSCycles(ii,:) = [Q(i) Gamma(j) iter];
        SteadyStateDensityTemp(ii,:) = y(end,:);
        

    end
    toc
    delete(poolobj);

%% Restructure results into 2D array
for ii = 1:length(Q)*length(Gamma)

        [j,i]=ind2sub([length(Gamma),length(Q)],ii);
        SteadyStateDensity(i,j,:) = SteadyStateDensityTemp(ii,:);
end

%% Save workspace

filename = sprintf("CyclePeriod=%.1f,S0=%1.e,V0=%1.e,p_L=%.1f,p_V=%.1f.mat",CyclePeriod,S0,V01,p_L,p_V);

save(filename);

end

