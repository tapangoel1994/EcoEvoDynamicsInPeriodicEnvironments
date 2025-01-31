%% This script generates data for Figures 4 and 5 for the revised manuscript.
%For a range of cycle periods and filtration conditions, this script finds
%the steady state viral genome densities across different viral strategies
%and for a fixed gamma=gamma*, it calculates the PIPs for those cycle
%period and filtration conditions.

%% Note: This code uses a parfor loop so you need the matlab parallel computing toolbox to run this script

%% Date last modified: 1/31/2025
%% Author: Tapan Goel

clear all;
close all;

%% load fixed life history parameters and colorscheme
addpath('utils/','lib/','../Data/');
colorpalette;
fixedparameters_highdeathrate;

%% Define parameters
Gamma = 10.^(-6:.1:-1);
P = 0:.02:1;

NumNodes = 64; %number of nodes used for parallel computing (if needed)

CyclePeriodList = 24;%[4, 8, 12, 16, 20, 24, 36, 48, 96];

q_LV = [.1 .1];%[0.025 0;0.05 0;0.1 0;0.2 0;0.4 0;...
       % 0 0.05;0 0.1;0 0.2];

    
%% Load steady state density matrices (either from a datafile if it already exists or by running the utils\PopulationSteadyStateFunction if the data doesnt exist
for cycleperiodindex = 1:length(CyclePeriodList)
    for filterindex = 1:size(q_LV,1)
        filename = sprintf('../Data/SteadyState_CyclePeriod=%d,d=%.2f,q_L=%.3f,q_V=%.3f.mat',CyclePeriodList(cycleperiodindex), ...
            params.d_S,q_LV(filterindex,1),q_LV(filterindex,2));
        if isfile(filename)
            load(filename,"SteadyStateDensity","SSCycles");
            SteadyState{cycleperiodindex,filterindex} = SteadyStateDensity;
            CyclesSteadyState{cycleperiodindex,filterindex} = SSCycles;
        else
            Save.Flag = 1;
            Save.FileName = filename;
            [SteadyState{cycleperiodindex,filterindex}, CyclesSteadyState{cycleperiodindex,filterindex}] = PopulationSteadyStateFunction(CyclePeriodList(cycleperiodindex), ...
                q_LV(filterindex,1),q_LV(filterindex,2),MaxCycles,Gamma,P,NumNodes,Save,params);
        end
    end
end

%% Load or generate invasion matrices
OptimalStrategy = cell(length(CyclePeriodList),size(q_LV,1));
for cycleperiodindex = 1:length(CyclePeriodList)
    for filterindex = 1:size(q_LV,1)
        
        % Find optimal viral strategy
        SteadyState_temp = SteadyState{cycleperiodindex,filterindex};
        SteadyState_temp = squeeze(sum(SteadyState_temp(:,:,3:10),3));
        
        [r,c] = find(SteadyState_temp == max(SteadyState_temp,[],"all","linear"));
        j = max(r);
        i = max(c);
        OptimalStrategy{cycleperiodindex,filterindex} = [P(j) Gamma(i) max(SteadyState_temp,[],"all","linear")];
        
        if OptimalStrategy{cycleperiodindex,filterindex}(3) > 10*criticaldensitythreshold % doing invasion analysis for the optimal strategy only makes sense if the strategy produces more viral genomes than the lower threshold
            
            InvasionVariable = [P' Gamma(i)*ones(size(P'))];

            filename = sprintf('../Data/Invasion_CyclePeriod=%d,d=%.2f,Gamma=%.3e,q_L=%.3f,q_V=%.3f.mat',CyclePeriodList(cycleperiodindex), ...
                params.d_S,InvasionVariable(1,2),q_LV(filterindex,1),q_LV(filterindex,2));

            if isfile(filename)
                
                load(filename,"InvasionDensity", "InvasionMatrix", "CyclesToInvasion");
                Invasion{cycleperiodindex,filterindex} = InvasionDensity;
                CyclesInvasion{cycleperiodindex,filterindex} = CyclesToInvasion;
                InvasionSuccessMatrix{cycleperiodindex,filterindex} = InvasionMatrix;
                
            else
                Save.Flag = 1;
                Save.FileName = filename;
                
                [Invasion{cycleperiodindex,filterindex}, CyclesInvasion{cycleperiodindex,filterindex}, InvasionSuccessMatrix{cycleperiodindex,filterindex}] = ...
                InvasionDynamics(CyclePeriodList(cycleperiodindex),q_LV(filterindex,1),q_LV(filterindex,2),InvasionVariable,MaxCycles,NumNodes,Save, params);
            end

        end

    end
end

