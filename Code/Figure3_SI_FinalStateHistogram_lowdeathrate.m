%% SI for figure 3 to compare total steady state viral genome densities across filtration conditions and strategies
%% Generates figure S4 

%% Date created:1/18/2025
%% Date modified: 1/31/2025
%% Author: Tapan Goel

close all; 
clear all;

%% load fixed parameters and color scheme
addpath('utils/');
addpath('lib/');
    colorpalette;
    fixedparameters_lowdeathrate;


%% Experimental parameters
Q_LV = [0 .2;.1 .1;.2 0];
params.T = 24;
params.t_vals = transpose(0:params.dt:params.T);

%% Viral strategies
ViralStrategies = [0 0;.92 .006;1 0]; %% each row contains the p-value and gamma-value for a viral straetgy
StrategyLabels = {'Obligately lytic','Temperate','Obligately lysogenic'};
FiltrateLabels = {'virions only','virions and lysogens','lysogens only'};

%% initial conditions
R0 = 100; %initial resource amount in ug/mL ( 500 mL flask)
S0 = 1e7; %Initial concentration of susceptibles in flask (per mL)
Va_0= 1e4; %initial concentration of virus in flask (per mL)
Vb_0 = 0;
x0 = [R0 S0 zeros(1,6) Va_0 Vb_0];

%% Calculate steady state for each strategy for each filter type
SteadyStates = zeros(size(Q_LV,1),size(ViralStrategies,1),length(x0));


for filterindex = 1:size(Q_LV,1)
    for strategyindex = 1:size(ViralStrategies,1)
        
        %% Set filtration condition
        q_L = Q_LV(filterindex,1);
        q_V = Q_LV(filterindex,2);
        TransferMatrix = diag([q_R q_S q_E q_E q_I q_I q_L q_L q_V q_V]);
        
        %% Set viral strategy
        params.p = [ViralStrategies(strategyindex,1) 0];
        params.gamma = [ViralStrategies(strategyindex,2) 0];

        % first cycle
        x0 = [R0 S0 zeros(1,6) Va_0 Vb_0]; %initial condition
        [t_vals, y] = ode113(@ODE_RSEILV_2Species, params.t_vals, x0, options, params); %growth cycle
        T = t_vals;
        Y = y;
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
            
            [t_vals, y] = ode113(@ODE_RSEILV_2Species, params.t_vals, x0, options, params); %growth cycle
             T = [T; params.T*(iter-1) + t_vals];
             Y = [Y;y];
            x0 = [R0 S0 zeros(1,8)] + y(end,:)*TransferMatrix; %transfer
        end
            
        %% Store steady state final condition
        SteadyStates(filterindex,strategyindex,:) = x0;
        
    end

end



SteadyStateViralGenomeDensities = sum(SteadyStates(:,:,3:10),3);

% Create a grouped bar chart
h = figure('Position',[767,599,793,639]);
bar(SteadyStateViralGenomeDensities', 'grouped');  % Transpose A to group bars by environment

% Labeling the plot
title('Total viral genome densities for different viral strategies','FontSize',16);
xlabel('Viral strategy','Interpreter','latex','FontSize',14);
ylabel('Total viral genome density (mL$^{-1}$)','Interpreter','latex','FontSize',14);
xticks(1:size(StrategyLabels,2));  % X-axis labels corresponding to the environments
xticklabels(StrategyLabels);  % Customize if needed

set(gca,'YScale','log')
set(gca,'YLim',[1e-4 1e10])
set(gca,'YMinorTick','off')
set(gca,'FontSize',14)

hold on;
yline(criticaldensitythreshold,':r');
% Add legend to show which bar corresponds to which strategy
lgd = legend(FiltrateLabels, 'Location', 'northeast','Interpreter','latex','FontSize',14);
lgd.Title.String = 'Filtrate:';       
            
%% Save Figure;
saveas(h,'../RevisedFigures/FigureS4.svg');       
       


       
             