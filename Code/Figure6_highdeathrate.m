% This script generates figure S13 of the manuscript

%% Note: The script generates a figure with annotations. The annotation placement has been modified in inkscape for better clarity.

%% Author: Tapan Goel
%% Last Modified: 1/31/2025

clear all;
close all;


addpath('lib/','utils/','../Data/');
VariableMaps;

run('fixedparameters_highdeathrate');
run('colorpalette');
StrategyNames = {'Obligately lytic','Temperate','Obligately lysogenic'};
%% Set stochastic realization parameters;
NumCycles = 100;
params.T = 24;
params.t_vals = 0:params.dt:params.T;
rng(1);

q_L = 0.1*ones(NumCycles,1);

q_V = 10.^(-1-5*rand(NumCycles,1)); %log uniform distribution over [1e-5, 1e-1]

%%Create figure;

f = figure('Position',[360,96,757,1187]);
t = tiledlayout(4,1);

nexttile(1);

plot(1:NumCycles,q_V,'LineStyle','-','LineWidth',1,'Color','k');
hold on;
scatter(1:5:NumCycles,q_V(1:5:NumCycles),15,'MarkerFaceColor','k','MarkerEdgeColor','k');
ylabel('$q_V$','FontSize',14,'FontWeight','bold','Interpreter','latex','Rotation',0);
set(gca,'YScale','log','YTick',10.^(-5:1:-1),'XTick',0:10:NumCycles,'box','off','YMinorTick','off','FontSize',12);
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
    nexttile(strategyindex+1);
       
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

nexttile(2);
    lgd = legend(StateVariables1D(2:end),'FontSize',12,'Interpreter','latex','Box','off','Location','best');

nexttile(4);
    xlabel('Number of cycles completed','FontSize',14,'FontWeight','bold','Interpreter','latex');
    



%% Add panel labels

panellabels = 'ABCD';

label_to_tilemap = [1 1;2 2;3 3;4 4];

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
saveas(f,'../RevisedFigures/FigureS13.svg');










