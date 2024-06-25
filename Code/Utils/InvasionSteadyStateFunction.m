%%Code takes one-host two-virus system and sweeps over q for a fixed gamma or gamma for a fixed q keeping
%other parameters fixed and does invasion analysis. Uses an RSEILV model with no MOI dependence. 
% Obtain appropriate steady state densities for
% every resident mutant pair and generate a heatmap. 

%%Date Created: 1/23/2024
%%Author: Tapan Goel

%% Inputs:
% CyclePeriod - Duration of Single Cycle in hours.
% p_L - fraction of lysogens being passaged between cycles.
% p_V - fraction of virions being passaged between cycles.
% InvasionVariable - a Nx2 array where each row is a (q,gamma) pair that
%                    defines a species.
% numNodes - number of nodes being used for parallel computing. Do not
%            exceed 32 (depending on cluster being used).
% SaveFlag - set to 1 if you want to save the output workspace to a
%           .matfile
% params as the varargin input - simulation parameters including life history traits and
% simulation parameters. if params is empty, the function assigns default
% values to the parameters internally.

%% Output:
% Script generates a .mat file named:
% "InvasionCyclePeriod=<CyclePeriod>,S0=<InitialHostDensity>,V0=<InitialViralDensity>,p_L=<p_L>,p_V=<p_V>.mat"
% that contains all the variables in the workspace from the simulation.
% Note that the output file has variables for the steady state values of
% each cell type for all strategies evaluated but does not have the
% dynamics for the strategies.

function [InvasionSteadyStateDensity, InvasionSSCycles] = InvasionSteadyStateFunction(CyclePeriod,p_L,p_V,InvasionVariable,numNodes,SaveFlag, varargin)

%% If life history and simulation parameters are not added as a function input, create parameter values
if nargin == 5

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

else 
    params = varargin{1}; %% if life history and simulation parameters were added as function input, assign them to the params variable.
end


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

InvasionSteadyStateDensity = zeros(length(InvasionVariable),length(InvasionVariable),10);
InvasionSSCycles = zeros(length(InvasionVariable),length(InvasionVariable));

%% initial conditions
R0 = 1e2; %initial resource amount in ug/mL ( 500 mL flask)
S0 = 1e7; %Initial concentration of susceptibles in flask (per mL)
V01= 1e4; %initial concentration of virus in flask (per mL)
V02 = 1e2;

%% Initiate parallel pool
poolobj = parpool(numNodes);

tic


parfor resident = 1:length(InvasionVariable)

        Params = params;
        Params.q = [InvasionVariable(resident,1) 0];
        Params.gamma = [InvasionVariable(resident,2) 0];
        %% Resident dynamics to steady state
        
        %First cycle
        x0 = [R0 S0 zeros(1,6) V01 V02];
        [t_vals, y] = ode113(@ODE_RSEILV_2Species, Params.t_vals, x0, options, Params);        
        % Cycles till steady state or till MaxCycles
        iter = 1;
        steadyrep = 0;
        x0 = [R0 S0 zeros(1,8)] + y(end,:)*TransferMatrix;
        % loop continues till:
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
        end
        
        % Add mutant and do invasions
        InvasionIC = x0 + [zeros(1,9) V02];
        InvasionforResident = zeros(length(InvasionVariable),10);
        InvasionCyclesforResident = zeros(length(InvasionVariable),1);
        for mutant = 1:length(InvasionVariable)
            if(mutant ~= resident)

                Params.q = [InvasionVariable(resident,1) InvasionVariable(mutant,1)];
                Params.gamma = [InvasionVariable(resident,2) InvasionVariable(mutant,2)];
                %First cycle
                [t_vals, y] = ode113(@ODE_RSEILV_2Species, Params.t_vals, InvasionIC, options, Params);        
                % Cycles till steady state or till MaxCycles
                iter = 1;
                steadyrep = 0;
                x0 = [R0 S0 zeros(1,8)] + y(end,:)*TransferMatrix;
                % loop continues till:
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
                end
            InvasionforResident(mutant,:) = x0;  
            InvasionCyclesforResident(mutant) = iter;
            
            end
        end
        InvasionSteadyStateDensity(resident,:,:) = InvasionforResident;
        InvasionSSCycles(resident,:) = InvasionCyclesforResident;

    end
    toc
    delete(poolobj);

%% Save workspace
if SaveFlag == 1
    filename = sprintf("InvasionCyclePeriod=%.1f,S0=%1.e,V0=%1.e,p_L=%.1f,p_V=%.1f.mat",CyclePeriod,S0,V01,p_L,p_V);
    save(filename);
end

end

