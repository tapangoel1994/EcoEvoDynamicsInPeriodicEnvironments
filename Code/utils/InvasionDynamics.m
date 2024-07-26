%% Code takes one-host two-virus system and sweeps over the invasion variable given by the user while keeping
% other parameters fixed and does invasion analysis. Uses an RSEILV model with no MOI dependence. 
% For a given resident(denoted by a)-mutant(denoted by b) pair, we first obtain the steady state densities for resident. Then, we add the mutant at
% threshold density and look at the first Invasion_Cycles = 10 cycles to see if the mutant frequency has
% increased or decreased to declare if the invasion is successful or not. If the answer is not conclusive after 10
% cycles, we run another 10 cycles.
% We repeat this procedure for every resident-mutant pair in the invasion
% variable vector.

%%Date Created: 07/16/2024
%%Author: Tapan Goel

%% Inputs:
% CyclePeriod - Duration of Single Cycle in hours.
% q_L - fraction of lysogens being passaged between cycles.
% q_V - fraction of virions being passaged between cycles.
% InvasionVariable - a Nx2 array where each row is a (p,gamma) pair that
%                    defines a species.
% numNodes - number of nodes being used for parallel computing. Do not
%            exceed 32 (depending on cluster being used).
% SaveFlag - set to 1 if you want to save the output workspace to a
%           .matfile
% params as the varargin input - simulation parameters including life history traits and
% simulation parameters. if params is empty, the function assigns default
% values to the parameters internally.

%% Output:
% InvasionDensity - 3D array of size
%                   length(InvasionVariable)xlength(InvasionVariable)x10.
%                   Y = InvasionDensity(i,j:) contains the density vector of the one host two virus system for resident with traits InvasionVariable(i,:) 
%                   and mutant with traits InvasionVariable(j,:), at the beginning of the cycle where invasion success/failure is determined.
%                   The vector itself is Y =
%                   [R,S,E_a,E_b,I_a,I_b,L_a,L_b,V_a,V_b].
% InvasionMatrix - 2D matrix of size length(InvasionVariable)xlength(InvasionVariable).
%                  InvasionMatrix(i,j) contains a value that determines
%                  whether the invasion by mutant with traits
%                  InvasionVariable(j,:) was successful in a system with
%                  resident with traits InvasionVariable(i,:).
%                  InvasionMatrix(i,j) = -2, if the resident is not present
%                  in the system at steady state to begin with (so the
%                  notion of invasion is not meaningful).
%                  InvasionMatrix(i,j) = -1, if the invasion fails.
%                  InvasionMatrix(i,j) = 0, if the resident and mutant have
%                  the same trait values (i.e., i = j).
%                  InvasionMatrix(i,j) = 1, if the invasion is successful
% CyclesToInvasion - length(InvasionVariable)xlength(InvasionVariable).
%                  CyclesToInvasion(i,j) is the number of cycles it takes 
%                  after the introduction of the mutant with traits
%                  InvasionVariable(j,:), in a system with
%                  resident with traits InvasionVariable(i,:), to determine
%                  whether the invasion was a success or a failure.
%
% This function also generates a .mat file named:
% "../Data/Invasion_CyclePeriod=<CyclePeriod>,S0=<InitialHostDensity>,V0=<InitialViralDensity>,q_L=<q_L>,q_V=<q_V>.mat"
% that contains all the variables in the workspace from the simulation.


%% Note: This function uses the parfor loop and therefore needs the MATLAB Parallel Computing Toolbox to run

function [InvasionDensity, InvasionMatrix, CyclesToInvasion] = InvasionDynamics(CyclePeriod,q_L,q_V,InvasionVariable,numNodes,SaveFlag, varargin)

%% If life history and simulation parameters are not added as a function input, create parameter values
if nargin == 6

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
    
    %% simulation parameters:
    criticaldensitythreshold = 1e-3; % concentration difference below which two concentrations are treated as identical in per mL
    params.flask_volume = 1/criticaldensitythreshold; %volume in mL
    params.dt = 1/30; % hours
    

else 
    params = varargin{1}; %% if life history and simulation parameters were added as function input, assign them to the params variable.
end


MaxCycles = 50000;
InvasionCycles = 10;
params.T = CyclePeriod; % hours
params.t_vals = transpose(0:params.dt:params.T); % time

%% filter parameters
q_R = 1;
q_S = 0;
q_E = 0;
q_I = 0;
q_L = q_L;
q_V = q_V;
q_Total = q_E+q_I+q_L+q_V;
TransferMatrix = diag([q_R q_S q_E q_E q_I q_I q_L q_L q_V q_V]);

%% Numerical method related parameters
options = odeset('AbsTol',1e-8,'RelTol',1e-8,'NonNegative',1:10); %Options for the ODE function call
criticaldensitythreshold = 1/params.flask_volume; % concentration difference below which two concentrations are treated as identical

InvasionDensity = zeros(length(InvasionVariable),length(InvasionVariable),10);
CyclesToInvasion = zeros(length(InvasionVariable),length(InvasionVariable));
InvasionMatrix = zeros(length(InvasionVariable),length(InvasionVariable));

%% initial conditions
R0 = 1e2; %initial resource amount in ug/mL ( 500 mL flask)
S0 = 1e7; %Initial concentration of susceptibles in flask (per mL)
Va_0= 1e4; %initial concentration of virus in flask (per mL)
Vb_0 = 10*criticaldensitythreshold; %Initial mutant concentration 10* threshold density

%% Initiate parallel pool
poolobj = parpool(numNodes);

tic


parfor resident = 1:length(InvasionVariable)
    
        %% Resident dynamics to steady state
        Params = params;
        Params.p = [InvasionVariable(resident,1) 0];
        Params.gamma = [InvasionVariable(resident,2) 0];

        %First cycle
        x0 = [R0 S0 zeros(1,6) Va_0 0];
        [t_vals, y] = ode113(@ODE_RSEILV_2Species, Params.t_vals, x0, options, Params);        
        % Cycles till steady state or till MaxCycles
        iter = 1;
        steadyrep = 0;
        x0 = [R0 S0 zeros(1,8)] + y(end,:)*TransferMatrix;
        % loop continues till:
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
        end
        
        %% Add mutant and do invasions
        InvasionIC = x0 + [0 0 0 q_E*Vb_0 0 q_I*Vb_0 0 q_L*Vb_0 0 q_V*Vb_0]./q_Total; %This condition adds mutant at a total density Vb_0 spread into the E,I,L,V classes according to their filtration ratios
        
        %define variables to circumvent scope rules for parfor
        InvasionforResidentDensity = zeros(length(InvasionVariable),10); %temporarily stores resident and mutant final state for a given resident
        InvasionMatrixforResident = zeros(length(InvasionVariable),1); %temporarily stores invasion success/failure status for every mutant for a fixed resident
        CyclesToInvasionforResident = zeros(length(InvasionVariable),1); %temporarily stores the number of cycles to invasion for every mutant given a fixed resident
    
        %loop over all possible mutants

        for mutant = 1:length(InvasionVariable)
            if(mutant ~= resident & sum(InvasionIC(3:2:end)) > 100*criticaldensitythreshold) 
                %invasion makes sense only if the mutant is different from the resident and if the resident is already present in the system at a high enough density

                Params.p = [InvasionVariable(resident,1) InvasionVariable(mutant,1)];
                Params.gamma = [InvasionVariable(resident,2) InvasionVariable(mutant,2)];
                %First cycle
                x0 = InvasionIC;
                InvasionTerminationFlag = 0; % set to 1 if there is clear success/failure determination about the invasion
                InvasionSeries = zeros(InvasionCycles+1,10);
                InvasionSeries(1,:) = x0;
                InvasionIterations = 0;

                while ~InvasionTerminationFlag %keep iterating till there is a clear success/failure determination about the invasion
                    
                    for ii = 1:InvasionCycles %run InvasionCycles = 10 growth cycles with the resident and the mutant in the system
                        [t_vals, y] = ode113(@ODE_RSEILV_2Species, Params.t_vals, x0, options, Params);
                        x0 = [R0 S0 zeros(1,8)] + y(end,:)*TransferMatrix;
                        InvasionSeries(ii+1,:) = x0;
                        InvasionIterations = InvasionIterations+1;
                    end
                    
                    if (sum(InvasionSeries(end,4:2:end),2) > Vb_0) & sum(InvasionSeries(end,4:2:end),2) > sum(InvasionSeries(end-1,4:2:end),2) %Invasion successful
                       InvasionTerminationFlag = 1;
                       InvasionMatrixforResident(mutant) = 1;
                    elseif (sum(InvasionSeries(end,4:2:end),2) < Vb_0) & sum(InvasionSeries(end,4:2:end),2) < sum(InvasionSeries(end-1,4:2:end),2) %Invasion failure     
                        InvasionTerminationFlag = 1;
                        InvasionMatrixforResident(mutant) = -1;
                    end
                end
                InvasionforResidentDensity(mutant,:) = InvasionSeries(end,:);  
                CyclesToInvasionforResident(mutant) = InvasionIterations;

            elseif (mutant ~= resident & sum(InvasionIC(3:2:end)) < 100*criticaldensitythreshold)
                InvasionforResidentDensity(mutant,:) = -1*ones(size(InvasionforResidentDensity(mutant,:)));
                CyclesToInvasionforResident(mutant) = 0;
                InvasionMatrixforResident(mutant) = -2;
            else
                InvasionMatrixforResident(mutant) = 0;
            end
        end
        InvasionDensity(resident,:,:) = InvasionforResidentDensity;
        CyclesToInvasion(resident,:) = CyclesToInvasionforResident;
        InvasionMatrix(resident,:) = InvasionMatrixforResident;
    end
    toc
    delete(poolobj);

    %% If a diagonal index falls inside the nonfeasible zone, add it to the non-feasible zone
        for diagindex = 2:length(InvasionVariable)-1
            if InvasionMatrix(diagindex,diagindex-1) == -2 & InvasionMatrix(diagindex,diagindex+1) == -2
                InvasionMatrix(diagindex,diagindex) = -2;
            end
        end
    
        if InvasionMatrix(1,2) == -2
            InvasionMatrix(1,1) = -2;
        end
        if InvasionMatrix(end,end-1) == -2
            InvasionMatrix(end,end) = -2;
        end
%% Save workspace
if SaveFlag == 1
    if ~isfolder('../Data/')
        mkdir('../Data/');
    end
    filename = sprintf("../Data/Invasion_CyclePeriod=%.1f,S0=%1.e,V0=%1.e,q_L=%.1f,q_V=%.1f.mat",CyclePeriod,S0,Va_0,q_L,q_V);
    save(filename);
end

end

