%%Code takes one-host one-virus system and sweeps over (p,gamma) keeping
%other parameters fixed. Uses an RSEILV model with no MOI dependence. 
% Obtain appropriate steady state densities for
% every (p,gamma) pair and generate a heatmap. 

%%Date Created: 1/23/2024
%%Author: Tapan Goel

%% Inputs:
% CyclePeriod - Duration of Single Cycle in hours.
% q_L - fraction of lysogens being passaged between cycles.
% q_V - fraction of virions being passaged between cycles.
% Gamma - array of Gamma values at which steady states are to be evaluated.
% P - array of P values at which steady states are to be evaluated. 
% numNodes - number of nodes being used for parallel computing. Do not
% exceed 32 (depending on cluster being used).
% SaveFlag - set to 1 if you want to save the output workspace to a
%           .matfile
% params as the varargin input - simulation parameters including life history traits and
% simulation parameters. if params is empty, the function assigns default
% values to the parameters internally.

%% Output:
% Script generates a .mat file named:
% "SteadyState_CyclePeriod=<CyclePeriod>,S0=<InitialHostDensity>,V0=<InitialViralDensity>,q_L=<q_L>,q_V=<q_V>.mat"
% that contains all the variables in the workspace from the simulation.
% Note that the output file has variables for the steady state values of
% each cell type for all strategies evaluated but does not have the
% dynamics for the strategies.

function [SteadyStateDensity, SSCycles] = PopulationSteadyStateFunction(CyclePeriod,q_L,q_V,Gamma,P,numNodes,SaveFlag,varargin)

%% If life history and simulation parameters are not added as a function input, create parameter values
if nargin == 7
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
    params.p = [0 0];
    params.gamma = [0 0];
    
    
    %% Simulation parameters:
    
    params.dt = 1/30; % hours
    MaxCycles = 50000; % Max number of cycles to steady state before while loop terminates
    InvasionCycles = 10;% Number of cycles in each set to evaluate transients during invasion
    criticaldensitythreshold = 1e-3; % concentration difference below which two concentrations are treated as identical in per mL
    params.flask_volume = 1/criticaldensitythreshold; %volume in mL      
else 
    params = varargin{1}; %% if life history and simulation parameters were added as function input, assign them to the params variable.
end

params.T = CyclePeriod; % hours
params.t_vals = transpose(0:params.dt:params.T); % time



%% filter parameters
q_R = 1;
q_S = 0;
q_E = 0;
q_I = 0;
q_L = q_L;
q_V = q_V;
TransferMatrix = diag([q_R q_S q_E q_E q_I q_I q_L q_L q_V q_V]);

%% Numerical method related parameters
options = odeset('AbsTol',1e-8,'RelTol',1e-8,'NonNegative',1:10); %Options for the ODE function call

MaxCycles = 50000; % Max number of cycles to steady state before while loop terminates
InvasionCycles = 10;% Number of cycles in each set to evaluate transients during invasion
criticaldensitythreshold = 1e-3; % concentration difference below which two concentrations are treated as identical in per mL
params.flask_volume = 1/criticaldensitythreshold; %volume in mL      


SteadyStateDensityTemp = zeros(length(P)*length(Gamma),10);
SteadyStateDensity = zeros(length(P),length(Gamma),10);
SSCyclesTemp = zeros(length(P)*length(Gamma),3); %% Variable to store p,gamma,#cyclestoSteadyState
SSCycles = zeros(length(P),length(Gamma));

%% initial conditions
R0 = 1e2; %initial resource amount in ug/mL ( 500 mL flask)
S0 = 1e7; %Initial concentration of susceptibles in flask (per mL)
Va_0= 1e4; %initial concentration of virus in flask (per mL)
Vb_0 = 0;

%% Initiate parallel pool
poolobj = parpool(numNodes);

tic
parfor ii = 1:length(P)*length(Gamma)
        Params = params;
        [j,i]=ind2sub([length(Gamma),length(P)],ii);
        Params.p = [P(i) 0];
        Params.gamma = [Gamma(j) 0];
        
        %% First cycle
        x0 = [R0 S0 zeros(1,6) Va_0 0];
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
        
        while ( (sum(abs(x0 - y(1,:)) > criticaldensitythreshold) >= 1) || steadyrep < 10) && iter < MaxCycles+10
            
            if((sum(abs(x0 - y(1,:)) > criticaldensitythreshold) < 1))
                        steadyrep = steadyrep+1;
            end
        
        iter = iter+1;
        
        [t_vals, y] = ode113(@ODE_RSEILV_2Species, Params.t_vals, x0, options, Params);
        x0 = [R0 S0 zeros(1,8)] + y(end,:)*TransferMatrix;
        %TimeSeries = [TimeSeries;y];
        %PlotTimeSeries_2species(t_vals,y);
        %drawnow;
        end
        SSCyclesTemp(ii,:) = [P(i) Gamma(j) iter];
        SteadyStateDensityTemp(ii,:) = x0;
        

    end
    toc
    delete(poolobj);

%% Restructure results into 2D array
for ii = 1:length(P)*length(Gamma)

        [j,i]=ind2sub([length(Gamma),length(P)],ii);
        SteadyStateDensity(i,j,:) = SteadyStateDensityTemp(ii,:);
        SSCycles(i,j) = SSCyclesTemp(ii,3);
end

%% Save workspace
if SaveFlag == 1
     if ~isfolder('..\Data\')
        mkdir('..\Data\');
    end
    filename = sprintf("..\\Data\\SteadyState_CyclePeriod=%.1f,S0=%1.e,V0=%1.e,q_L=%.1f,q_V=%.1f.mat",CyclePeriod,S0,Va_0,q_L,q_V);
    save(filename);
end

end

