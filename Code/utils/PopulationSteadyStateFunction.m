%% Code takes one-host one-virus system and sweeps over (p,gamma) keeping
% other parameters fixed. Uses an RSEILV model with no MOI dependence. 
% Obtain appropriate steady state densities for
% every (p,gamma) pair, store the steady state values and generate a heatmap. 

%%Date Created: 12/5/2024
%%Date updated: 1/9/2025
%%Author: Tapan Goel

%% Inputs:
% CyclePeriod - Duration of single growth cycle in hours.
% q_L - fraction of lysogens being passaged between cycles.
% q_V - fraction of virions being passaged between cycles.
% Gamma - array of Gamma values at which steady states are to be evaluated.
% P - array of P values at which steady states are to be evaluated. 
% numNodes - number of nodes being used for parallel computing. Do not
% exceed 32 (depending on cluster being used).
% Save - If Save.Flag == 1, save the output workspace to a matfile
% named Save.FileName.
% params as the varargin input - simulation parameters including life history traits and
% simulation parameters. if params is empty, the function assigns default
% values to the parameters internally.

%% Output:
% SteadyStateDensity - 3D array of size length(P)xlenth(Gamma)x10. 
%                      Y = SteadyStateDensity(i,j,:) contains the steady
%                       state initial viral density vector for P(i),Gamma(j).
%                       The vector itself is Y =
%                       [R,S,E_a,0,I_a,0,L_a,0,V_a,0]
% SSCycles - 2D matrix of size length(P)xlength(Gamma). SSCycles(i,j) is
% the number of cycles it takes for the system to reach steady state for
% P(i), Gamma(j).


% Note that the output matfile (if saved) has variables for the steady state values of
% each cell type for all strategies evaluated but does not have the
% dynamics for the strategies.

%% Note: This function uses the parfor loop and therefore needs the MATLAB Parallel Computing Toolbox to run

function [SteadyStateDensity, SSCycles] = PopulationSteadyStateFunction(CyclePeriod_Input,q_L,q_V,MaxCycles_Input,Gamma,P,numNodes,Save,varargin)

%% If life history and simulation parameters are not added as a function input, create parameter values
if isempty(varargin)
    %% Life history parameters (units of hours, micrograms and mL). 
    addpath('..\lib\')
    fixedparameters;
    params.p = [0 0];
    params.gamma = [0 0];
   
else 
    params = varargin{1}; %% if life history and simulation parameters were added as function input, assign them to the params variable.
end

%% Simulation parameters added as inputs to the function
params.T = CyclePeriod_Input; % hours
params.t_vals = transpose(0:params.dt:params.T); % time vector
MaxCycles = MaxCycles_Input;


%% filter parameters
q_R = 0;
q_S = 0;
q_E = 0;
q_I = 0;
q_L = q_L;
q_V = q_V;
TransferMatrix = diag([q_R q_S q_E q_E q_I q_I q_L q_L q_V q_V]);

%% Numerical method related parameters
options = odeset('AbsTol',1e-8,'RelTol',1e-8,'NonNegative',1:10); %Options for the ODE function call
criticaldensitythreshold = 1e-3; % concentration difference below which two concentrations are treated as identical in per mL
params.flask_volume = 1/criticaldensitythreshold; %volume in mL      


SteadyStateDensityTemp = zeros(length(P)*length(Gamma),10); %%stores steady state density vectors temporarily as a 2D matrix (to make it easier to use parfor)
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
parfor ii = 1:length(P)*length(Gamma)  %%parallel loop over all possible pairs of P and Gamma

        Params = params;
        [j,i]=ind2sub([length(Gamma),length(P)],ii); %%to figure out how the parallel pool index ii corresponds to the indices of P and Gamma
        Params.p = [P(i) 0];
        Params.gamma = [Gamma(j) 0];
        
        %% First cycle
        x0 = [R0 S0 zeros(1,6) Va_0 0]; %initial condition
        [t_vals, y] = ode113(@ODE_RSEILV_2Species, Params.t_vals, x0, options, Params); %growth cycle
                
        %% Cycles till steady state or till MaxCycles
        iter = 1;
        steadyrep = 0;
        
        x0 = [R0 S0 zeros(1,8)] + y(end,:)*TransferMatrix; %first filtration
        
        %% while loop continues till:
        %1. The initial population for the new epoch becomes (almost) equal to the
        %initial population for the previous epoch for atleast 10 cycles OR,
        %2. The total number of cycles hits the max number of cycles.
        
        while ( (sum(abs(x0 - y(1,:)) > criticaldensitythreshold) >= 1) || steadyrep < 10) && iter < MaxCycles+10 
            
            if((sum(abs(x0 - y(1,:)) > criticaldensitythreshold) < 1))
                        steadyrep = steadyrep+1;
            end
        
            iter = iter+1;
            
            [t_vals, y] = ode113(@ODE_RSEILV_2Species, Params.t_vals, x0, options, Params); %growth cycle
            x0 = [R0 S0 zeros(1,8)] + y(end,:)*TransferMatrix; %transfer
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
if Save.Flag == 1
    if ~isempty(Save.FileName)
        save(Save.FileName);
    else
        error('No filename provided');
    end
end

end
