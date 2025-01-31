%% Script runs each strategy through 100 realizations of 100 cycles each with different levels of fixed qL and different lowerbound levels of stochastic qV to generate heatmaps.


%% Author: Tapan Goel
%% Last Modified: 1/31/2025

clear all;
close all;

%% Load variable maps and fixed parameters
addpath('lib/','utils/','../Data/');
VariableMaps;

deathrate =  'highdeathrate';

run(['fixedparameters_' deathrate]);


%% Set ensamble parameters
QL_LowEnd = 1:.25:13;
QV_LowEnd = 1:.25:13;

NumSims = 100;
NumCycles = 100;
CyclePeriod = 24; %in hours

%% Set initial condition for virus
x0 = [R0 S0 zeros(1,6) Va_0 0]; 

poolobj = parpool(64);

parfor ii = 1:length(QL_LowEnd)*length(QV_LowEnd)*NumSims
    
    Params = params;
    [qLidx, qVidx, simnum] = ind2sub([length(QL_LowEnd), length(QV_LowEnd), NumSims], ii);
    %% Generate parameters to scan over;
    
    rng('shuffle','twister');
    
    q_L = -1 - (QL_LowEnd(qLidx)-1)*rand(NumCycles,1); % generate uniform random numbers in range [-QL_LowEnd(qLidx), -1];
    q_V = -1 - (QV_LowEnd(qVidx)-1)*rand(NumCycles,1); % generate uniform random numbers in range [-QV_LowEnd(qVidx), -1];

    q_L = 10.^q_L; %generate log uniform random numbers for q_L
    q_V = 10.^q_V; %generate log uniform random numbers for q_V

    CyclesToTerminationForStrategy = zeros(1,3);

    for strategyindex = 1:3
    
        Params.p = [PGamma.(Strategies{strategyindex})(1) 0];
        Params.gamma = [PGamma.(Strategies{strategyindex})(2) 0];
        Params.T = CyclePeriod;
        Params.t_vals = 0:Params.dt:Params.T;
        CyclesToTerminationForStrategy(strategyindex) = StochasticRealization(NumCycles,Params,criticaldensitythreshold,q_R, q_S,q_E,q_I,q_L,q_V,x0,options);
        
    end
    
    CyclesToTerminationTemp(ii,:) = CyclesToTerminationForStrategy;

end

delete(poolobj);

for ii = 1:length(QL_LowEnd)*length(QV_LowEnd)*NumSims

    [qLidx, qVidx, simnum] = ind2sub([length(QL_LowEnd), length(QV_LowEnd), NumSims], ii);

    CyclesToTermination(:,simnum,qVidx,qLidx) = CyclesToTerminationTemp(ii,:);

end

% Initialize the table with appropriate variable names
CyclesToTerminationTable = table('Size', [0, 5], ...
    'VariableTypes', {'double', 'double', 'double', 'string', 'double'}, ...
    'VariableNames', {'QL_LowEnd', 'QV_LowEnd', 'SimulationNumber', 'Strategy', 'CyclesToTermination'});

% Loop through each combination of QV_LowEnd, QL, and Simulation Number
for qLidx = 1:length(QL_LowEnd)
    for qVidx = 1:length(QV_LowEnd)
        for simnum = 1:NumSims
            for strategyindex = 1:3
                % Extract the CyclesToTermination value for the current combination
                cycles = CyclesToTermination(strategyindex, simnum, qVidx, qLidx);
                
                % Convert the strategy name to a string (if it's a character vector)
                strategyName = string(Strategies{strategyindex});
                
                % Create a new row for the table
                newRow = table(QL_LowEnd(qLidx), QV_LowEnd(qVidx), simnum, strategyName, cycles, ...
                    'VariableNames', {'QL_LowEnd', 'QV_LowEnd', 'SimulationNumber', 'Strategy', 'CyclesToTermination'});
                
                % Append the new row to the table
                CyclesToTerminationTable = [CyclesToTerminationTable; newRow];
            end
        end
    end
end


filename = ['Stochastic_QL_Stochastic_QV_' deathrate '.mat'];

save(filename);

