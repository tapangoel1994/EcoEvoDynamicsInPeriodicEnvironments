% This script generates figure 6 of the manuscript

%% Note: The script generates a figure with annotations. The annotation placement has been modified in inkscape for better clarity.
%% Note: Run script Stochastic_QL_Stochastic_QV_lowdeathrate.m & script Stochastic_QL_Stochastic_QV_highdeathrate.m before running this script.

%% Author: Tapan Goel
%% Last Modified: 1/31/2025

clear all;
close all;


addpath('lib/','utils/','../Data/');
VariableMaps;

run('fixedparameters_lowdeathrate');
run('colorpalette');
StrategyNames = {'Obligately lytic','Temperate','Obligately lysogenic'};
%% Set stochastic realization parameters;
NumCycles = 100;
params.T = 24;
params.t_vals = 0:params.dt:params.T;
rng(1);

q_L = 0.1*ones(NumCycles,1);

q_V = 10.^(-1-10*rand(NumCycles,1)); %log uniform distribution over [1e-10, 1e-1]

%%Create figure;

f = figure('Position',[360,96,1178,1187]);
t = tiledlayout(4,2);

nexttile(1);

plot(1:NumCycles,q_V,'LineStyle','-','LineWidth',1,'Color','k');
hold on;
scatter(1:5:NumCycles,q_V(1:5:NumCycles),15,'MarkerFaceColor','k','MarkerEdgeColor','k');
ylabel('$q_V$','FontSize',14,'FontWeight','bold','Interpreter','latex','Rotation',0);
set(gca,'YScale','log','YTick',10.^(-10:2:-1),'XTick',0:10:NumCycles,'box','off','YMinorTick','off','FontSize',12);
hold off;

for strategyindex = 1:length(Strategies)

    params.p = [PGamma.(Strategies{strategyindex})(1) 0];
    params.gamma = [PGamma.(Strategies{strategyindex})(2) 0];

    T = [];
    Y = [];
    

    x0 = [R0 S0 zeros(1,6) Va_0 0]; % intial condition
    CycleEndStates = x0;

    for cyclenum = 1:NumCycles

        [t,y] = ode113(@ODE_RSEILV_2Species, params.t_vals, x0, options, params); % run cycle
      
        T = [T; params.T*(cyclenum)-1 + t];%append cycle dynamics 
        Y = [Y;y];%append cycle dynamics
    
        TransferMatrix = diag([q_R q_S q_E q_E q_I q_I q_L(cyclenum) q_L(cyclenum) q_V(cyclenum) q_V(cyclenum)]); % update transfer matrix
        
        x0 = [R0 S0 zeros(1,8)] + y(end,:)*TransferMatrix; %update initial condition for next cycle.
        
        x0(x0  < criticaldensitythreshold) = 0; %if something is below threshold, set it to zero

        CycleEndStates(cyclenum+1,:) = y(end,:);

    end
    
    % plot cycle to cycle dynamics
    nexttile(2*strategyindex+1);
       
        for statevaridx = 2:6
            semilogy(0:NumCycles,CycleEndStates(:,StateVariableMap1D.(StateVariables1D{statevaridx})),'LineWidth',1,'Color',linecolors.(StateVariables1D{statevaridx}));
            hold on;               
        end
        
        for statevaridx = 2:6
            semilogy(0:5:NumCycles,CycleEndStates(1:5:end,StateVariableMap1D.(StateVariables1D{statevaridx})),'LineWidth',1,'Color',linecolors.(StateVariables1D{statevaridx}),'LineStyle','none', ...
                'Marker',marker.(Strategies{strategyindex}));
            hold on;               
        end

        yline(criticaldensitythreshold,':k','LineWidth',0.5);

        set(gca,'XLim',[0 NumCycles],'YLim', [1e-4 1e10],'YTick',10.^(-4:2:10),'XTick',0:10:NumCycles,'YMinorTick','off','Box','off','FontSize',12);

        ylabel('Density (mL$^{-1}$)','FontSize',14,'FontWeight','bold','Interpreter','latex');

        text(70,1e9,StrategyNames{strategyindex},'Interpreter','latex','FontSize',12);
    
end

nexttile(3);
    lgd = legend(StateVariables1D(2:end),'FontSize',12,'Interpreter','latex','Box','off','Location','best');

nexttile(7);
    xlabel('Number of cycles completed','FontSize',14,'FontWeight','bold','Interpreter','latex');
    


%% Make survival heatmaps

load('../Data/Stochastic_QL_Stochastic_QV_lowdeathrate.mat');

% Convert cycles to termination table to heatmap
CyclesToTerminationTable.QL_LowEnd = 10.^-CyclesToTerminationTable.QL_LowEnd;
CyclesToTerminationTable.QV_LowEnd = 10.^-CyclesToTerminationTable.QV_LowEnd;

QL_LowEnd = unique(CyclesToTerminationTable.QL_LowEnd);
QV_LowEnd = unique(CyclesToTerminationTable.QV_LowEnd);

for strategyidx = 1:length(Strategies)

    currData = CyclesToTerminationTable(strcmp(CyclesToTerminationTable.Strategy,Strategies{strategyidx}),:);
        
    for i = 1:length(QL_LowEnd)
        for j = 1:length(QV_LowEnd)
            

            SurvivalProbability_Matrix(i,j,strategyidx) = ...
                height(currData(currData.QL_LowEnd == QL_LowEnd(i) & currData.QV_LowEnd == QV_LowEnd(j) & currData.CyclesToTermination == NumCycles,:)); %% gives number of realizations where NCycles are reached

        end
    end
end

SurvivalProbability_Matrix = SurvivalProbability_Matrix/NumCycles; %% Convert counts to probability


nexttile(2,[2,1]);

  imagesc(QL_LowEnd,QV_LowEnd,squeeze(SurvivalProbability_Matrix(:,:,2))');
    xlabel('Lower limit of $q_L$','Interpreter','latex','FontSize',14);
    ylabel('Lower limit of $q_V$','Interpreter','latex','FontSize',14);
    set(gca,'PlotBoxAspectRatio',[1 1 1]);
    set(gca,'XScale','log','YScale','log','CLim',[0 1],'XTick',10.^(-13:2:-1),'YTick',10.^(-13:2:-1),'FontSize',12,'TickLabelInterpreter','latex');
    colormap("parula");
    c = colorbar;
    c.Label.String = 'Probability of survival';
    c.Label.Rotation = -90;
    c.TickLabelInterpreter = 'latex';
    c.Label.Interpreter = 'latex';
    axis xy;

    hold on;

    c1 = contourc(QL_LowEnd, QV_LowEnd, squeeze(SurvivalProbability_Matrix(:, :, 1))',[.05 .05]);
    c2 = contourc(QL_LowEnd, QV_LowEnd, squeeze(SurvivalProbability_Matrix(:, :, 1))',[.95 .95]);
    
    c1 = c1(:,2:end);
    c2 = c2(:,2:end);

    y_low = mean(c1(2,c1(2,1:end) <= max(QV_LowEnd))); % have to do this additional filtering of c1 because sometimes contourc C generates a garbage datapoint to close the contour.
    y_high = mean(c2(2,c2(2,1:end) <= max(QV_LowEnd)));

    yline(y_low,':','LineWidth',2.5,'Color',[1 1 1]);
    yline(y_high,'--','LineWidth',2.5,'Color',[1 1 1]);

    title('Cell death rate = 0.04 hr$^{-1}$','Interpreter','latex','FontSize',14);


load('Stochastic_QL_Stochastic_QV_highdeathrate.mat');

% Convert cycles to termination table to heatmap
CyclesToTerminationTable.QL_LowEnd = 10.^-CyclesToTerminationTable.QL_LowEnd;
CyclesToTerminationTable.QV_LowEnd = 10.^-CyclesToTerminationTable.QV_LowEnd;

QL_LowEnd = unique(CyclesToTerminationTable.QL_LowEnd);
QV_LowEnd = unique(CyclesToTerminationTable.QV_LowEnd);

for strategyidx = 1:length(Strategies)

    currData = CyclesToTerminationTable(strcmp(CyclesToTerminationTable.Strategy,Strategies{strategyidx}),:);
        
    for i = 1:length(QL_LowEnd)
        for j = 1:length(QV_LowEnd)
            

            SurvivalProbability_Matrix(i,j,strategyidx) = ...
                height(currData(currData.QL_LowEnd == QL_LowEnd(i) & currData.QV_LowEnd == QV_LowEnd(j) & currData.CyclesToTermination == NumCycles,:)); %% gives number of realizations where NCycles are reached

        end
    end
end

SurvivalProbability_Matrix = SurvivalProbability_Matrix/NumCycles; %% Convert counts to probability


nexttile(6,[2,1]);

  imagesc(QL_LowEnd,QV_LowEnd,squeeze(SurvivalProbability_Matrix(:,:,2))');
    xlabel('Lower limit of $q_L$','Interpreter','latex','FontSize',14);
    ylabel('Lower limit of $q_V$','Interpreter','latex','FontSize',14);
    
    set(gca,'PlotBoxAspectRatio',[1 1 1]);
    set(gca,'XScale','log','YScale','log','CLim',[0 1],'XTick',10.^(-13:2:-1),'YTick',10.^(-13:2:-1),'FontSize',12,'TickLabelInterpreter','latex');
    c = colorbar;
    c.Label.String = 'Probability of survival';
    c.Label.Rotation = -90;
    c.TickLabelInterpreter = 'latex';
    c.Label.Interpreter = 'latex';
    axis xy;

    hold on;

   c1 = contourc(QL_LowEnd, QV_LowEnd, squeeze(SurvivalProbability_Matrix(:, :, 1))',[.05 .05]);
    c2 = contourc(QL_LowEnd, QV_LowEnd, squeeze(SurvivalProbability_Matrix(:, :, 1))',[.95 .95]);
    
    c1 = c1(:,2:end);
    c2 = c2(:,2:end);

    y_low = mean(c1(2,c1(2,1:end) <= max(QV_LowEnd))); % have to do this additional filtering of c1 because sometimes contourc C generates a garbage datapoint to close the contour.
    y_high = mean(c2(2,c2(2,1:end) <= max(QV_LowEnd)));

    yline(y_low,':','LineWidth',2.5,'Color',[1 1 1]);
    yline(y_high,'--','LineWidth',2.5,'Color',[1 1 1]);

    title('Cell death rate = 0.2 hr$^{-1}$','Interpreter','latex','FontSize',14);

%% Add panel labels

panellabels = 'ABCDEF';

label_to_tilemap = [1 1;2 3;3 5;4 7;5 2;6 6];

for index = 1:length(panellabels)

    if index == 5 | index == 6

        nexttile(label_to_tilemap(index,2),[2 1]);
    else
        nexttile(label_to_tilemap(index,2))
    end

    xloc = get(gca,'XTick');

    xloc = min(xloc);

    yloc = get(gca,'YTick');

    yloc = max(yloc);

    text(xloc,yloc, sprintf("(%s)",panellabels(index)),'FontSize',14,'FontWeight','bold');

end


%% Save Figure
saveas(f,'../RevisedFigures/Figure6.svg');









