%% Script generates figures S1-S3 and S6-8
%% Date created: 6/13/2024
%% Date modified: 1/31/2025
%% Author: Tapan Goel

close all; 
clear all;

%% Figure making conditions
SaveFigureFlag = 0; % if SaveFigureFlag is 1, save the figure as eps and as png.
PlotAllDeathRates = 1; %if this variable is 1, generate figures for both low and high death rates, otherwise only generate figure for low death rate.

%% Add file paths
addpath('utils/','lib/','../Data/');
colorpalette;
VariableMaps;

%% Set Cycle period and filtration conditions
Q_LV = [0 .2;.1 .1;.2 0];
CyclePeriod = 24;
NCycles = 40;
%% Set strategies and labels

StrategyLabels = {'Obligately lytic','Temperate','Obligately lysogenic'};
ParameterFiles = {'lowdeathrate','highdeathrate'};
Filtrate = {'virions only','virions and lysogens','lysogens only'};

%% Create maps:
 strategytotile = [1 3;6 8;11 13]; %strategytotile(strategyindex,1) maps to the tiles corresponding to the plot for that particular viral strategy
 statevariablemap = struct('R',1,'S',2,'E',3,'I',5,'L',7,'V',9); % maps state variables to indices of the state vector Y.
 statevariables = {'R','S','E','I','L','V'};
 plotlabel = ['A' 'B';'C' 'D';'E' 'F'];
 
for filtrationindex = 1:size(Q_LV,1)
    for celldeathrateindex = 1:(PlotAllDeathRates+1)

       h_main = figure('Position',[100 0 1400 800]);
       t_main = tiledlayout(3,5,'TileSpacing','loose','Padding','loose');
        %% Load parametervalues
        run(['fixedparameters_' ParameterFiles{celldeathrateindex}]);

        %% Simulation parameters
        params.T = CyclePeriod;
        params.t_vals = transpose(0:params.dt:params.T);

        %% initial conditions
        R0 = 100; %initial resource amount in ug/mL ( 500 mL flask)
        S0 = 1e7; %Initial concentration of susceptibles in flask (per mL)
        Va_0= 1e4; %initial concentration of virus in flask (per mL)
        Vb_0 = 0;

        %% filtration conditions
        q_L = Q_LV(filtrationindex,1);
        q_V = Q_LV(filtrationindex,2);
            
        TransferMatrix = diag([q_R q_S q_E q_E q_I q_I q_L q_L q_V q_V]);
    
        %% Generate dynamics for each strategy
        for strategyindex = 1:length(Strategies)

            params.p = [PGamma.(Strategies{strategyindex})(1) 0];
            params.gamma = [PGamma.(Strategies{strategyindex})(2) 0];
            
            T = [];
            Y = [];
            x0 = [R0 S0 zeros(1,6) Va_0 Vb_0];
            
            for iter = 1:NCycles
                [t_vals, y] = ode113(@ODE_RSEILV_2Species, params.t_vals, x0, options, params);
                T = [T; params.T*(iter-1) + t_vals];
                Y = [Y;y];
                x0 = [R0 S0 zeros(1,8)] + y(end,:)*TransferMatrix;
            end
            
            %% Plot first 3 cycles
            nexttile(t_main,strategytotile(strategyindex,1),[1 2]);
            
            yline(criticaldensitythreshold,':r','LineWidth',1); hold on;
            for statevariableindex = 1:length(statevariables)
                semilogy(T(1:3*length(t_vals)),Y(1:3*length(t_vals),statevariablemap.(statevariables{statevariableindex})), ... 
                    'LineWidth',2,'Color',linecolors.(statevariables{statevariableindex}),'LineStyle','-'); 
            end

            for statevariableindex = 1:length(statevariables)
                semilogy(T(1:length(t_vals):3*length(t_vals)+1),Y(1:length(t_vals):3*length(t_vals)+1,statevariablemap.(statevariables{statevariableindex})), ... 
                    'LineWidth',1.5,'Color',linecolors.(statevariables{statevariableindex}),'LineStyle','none','marker',marker.(Strategies{strategyindex}));
            end
                        
            xline(.99*params.T,'--','LineWidth',1.5,'Color',.8*[1 1 1]);
            xline(1.99*params.T,'--','LineWidth',1.5,'Color',.5*[1 1 1]);
            xline(2.99*params.T,'--','LineWidth',1.5,'Color',.2*[1 1 1]);
            ylim([1e-4 1e10]);
            xlim([0 3.1*params.T]);
            set(gca,'YScale','log','YMinorTick','off','Box','off','XTick',0:params.T/3:3*params.T,'XTickLabel',[],'YTick',10.^(-4:2:10));
            set(gca,'FontSize',14,'FontWeight','bold','FontName','Times');
            text(-13*params.T/24,5e10,sprintf("(%s)",plotlabel(strategyindex,1)),'FontSize',16,'FontWeight','bold');
            hold off;
            
            %% Plot cycle-to-cycle  dynamics for the beginning of a cycle
            nexttile(t_main,strategytotile(strategyindex,2),[1 3]);
            
            CycleTimePoints = T(1:length(t_vals):end);
            CycleStartDensities = Y(1:length(t_vals):end,:);
            
            for statevaridx = 1:6   %% for the purposes of plotting, set all variables that are zero, equal to critical densitythreshold/10 for the first instance of zero and set them to zero after
                locs = find(CycleStartDensities(:,statevariablemap.(statevariables{statevaridx})) <= criticaldensitythreshold,1);
                if ~isempty(locs) & locs < size(CycleStartDensities,1)
                    if all(CycleStartDensities(locs+1:end,statevariablemap.(statevariables{statevaridx})) < criticaldensitythreshold) | locs == 1
                        CycleStartDensities(locs,statevariablemap.(statevariables{statevaridx})) = criticaldensitythreshold/10;
                    end
                end
            end

            yline(criticaldensitythreshold,':r','LineWidth',1);hold on;

            for statevariableindex = 1:length(statevariables)
                semilogy(CycleTimePoints,CycleStartDensities(:,statevariablemap.(statevariables{statevariableindex})), ... 
                    'LineWidth',1.5,'Color',linecolors.(statevariables{statevariableindex}),'LineStyle','-');
            end
            
            for statevariableindex = 1:length(statevariables)
                semilogy(CycleTimePoints([1:4 6:2:end],1),CycleStartDensities([1:4 6:2:end],statevariablemap.(statevariables{statevariableindex})), ...
                    'Marker',marker.(Strategies{strategyindex}),'Color',linecolors.(statevariables{statevariableindex}),'LineStyle','none');
            end

            xline(.99*params.T,'--','LineWidth',1.5,'Color',.8*[1 1 1]);
            xline(1.99*params.T,'--','LineWidth',1.5,'Color',.5*[1 1 1]);
            xline(2.99*params.T,'--','LineWidth',1.5,'Color',.2*[1 1 1]);
            ylim([1e-4 1e10]);
            xlim([0 T(end)+0.1*params.T]);
            
            set(gca,'YScale','log','YMinorTick','off','Box','off','XTick',params.T*[0:1:3 5:2:NCycles],'XTickLabel',[],'YTick',10.^(-4:2:10),'YTickLabel',[]);
            set(gca,'FontSize',14,'FontWeight','bold','FontName','Times');
            text(-60*params.T/24,5e10,sprintf("(%s)",plotlabel(strategyindex,2)),'FontSize',16,'FontWeight','bold');
            hold off;

        end
        
        %% Add Tick and axis labels and title.

        nexttile(t_main,strategytotile(end,1),[1 2]);
        set(gca,'XTickLabel',0:params.T/3:3*params.T);
        xlabel('Time (hr)','FontSize',16,'FontWeight','bold','FontName','Times','Interpreter','latex');

        nexttile(t_main,strategytotile(end,2),[1 3]);
        set(gca,'XTickLabel',[0:1:3 5:2:NCycles]);
        xlabel('Number of cycles completed','FontSize',16,'FontWeight','bold','FontName','Times','Interpreter','latex');

        title(t_main,sprintf("Filtrate: %s",Filtrate{filtrationindex}),'FontSize',16,'FontWeight','bold','FontName','Times','Interpreter','latex');
        ylabel(t_main,'Density (mL$^{-1}$)','FontSize',16,'FontWeight','bold','FontName','Times','Interpreter','latex');

        %% Add legends
        nexttile(t_main,strategytotile(1,1),[1 2]);
        lgd = legend('',statevariables{:},'Box','off','Location','best');
        
        nexttile(t_main,strategytotile(1,2),[1 3]);
        hold on;
        lh(1) = plot(nan,nan,'LineStyle','-','Marker',marker.lytic,'LineWidth',1,'Color','k');
        lh(2) = plot(nan,nan,'LineStyle','-','Marker',marker.temperate,'LineWidth',1,'Color','k');
        lh(3) = plot(nan,nan,'LineStyle','-','Marker',marker.lysogenic,'LineWidth',1,'Color','k');
        lgd = legend(lh,Strategies,'Box','off','Location','best');


        %% Save figure

        if SaveFigureFlag == 1

            filename = sprintf('../RevisedFigures/Figure3_SI_%s_%s.eps',ParameterFiles{celldeathrateindex},Filtrate{filtrationindex});          
            saveas(h_main,[filename(1:end-4) '.svg']);%save as svg

        end
                
    end
end
