%% Code to simulate and plot single cycle for 3 different viral strategies over 24
%% hours, for two different cell death rates. This script generates Figure 2 in the manuscript.

%% Date created: 12/15/2024
%% Date updated: 1/9/2025
%% Author: Tapan Goel

close all; 
clear all;

%% load fixed parameters and color scheme
addpath('utils/','lib/');
fixedparameters_lowdeathrate; %run script to load parameters for the low cell death rate scenario
colorpalette;  %run script to load line colors, markers etc. 
VariableMaps; %maps variable names to array indices for population density vectors

SaveFigureFlag = 1; % set to 1 to save the figures generated

%% simulation parameters:
params.T = 24; % hours
params.t_vals = transpose(0:params.dt:params.T); % time


%% initial conditions
R0 = 100; %initial resource amount in ug/mL
S0 = 1e7; %Initial concentration of susceptibles in flask (per mL)
V01= 1e4; %initial concentration of virus in flask (per mL)
V02 = 0;
x0 = [R0 S0 zeros(1,6) V01 V02];

PGamma = {[0 0;.92 .006;1 0],[0 0;.5 .083;1 0]}; %% each row is a (p,gamma) for a particular strategy. The temperate strategy is the optimal strategy for 24hr cycle, q_L = .2, q_V = 0 filter.
%the first cell is the strategy for low death rate. the second cell array
%is strategy for high death rate
StrategyLabels = {'Obligately lytic','Temperate','Obligately lysogenic'};

celldeathrates = [0.04 0.2]; % per capita death rates of susceptible, exposed, lytically infected and lysogenic cells in /hr

h = figure('Position',1.0e+03 *[1.0000    0.0770    0.8883    1.1873]);

tiles = tiledlayout(size(PGamma{1},1)+1, length(celldeathrates),'TileSpacing','compact','Padding','compact');

for deathrateindex = 1:length(celldeathrates)

    for strategyindex = 1:size(PGamma{1},1)
        
    
        params.d_S = celldeathrates(deathrateindex);
        params.d_E = celldeathrates(deathrateindex);
        params.d_I = celldeathrates(deathrateindex);
        params.d_L = celldeathrates(deathrateindex);
        
        
            params.p = [PGamma{deathrateindex}(strategyindex,1) 0];
            params.gamma = [PGamma{deathrateindex}(strategyindex,2) 0];
        
            [t_vals,y] = ode113(@ODE_RSEILV_2Species,params.t_vals,x0,options,params); %simulate one growth cycle
            
            %plot the timeseries for the 1 growth cycle.
            nexttile(deathrateindex+2*(strategyindex-1));
            yline(criticaldensitythreshold,':r','LineWidth',1);
            hold on;
            semilogy(t_vals,y(:,1),'LineWidth',3,'Color',linecolors.R);  %Plot R
            semilogy(t_vals,y(:,2),'LineWidth',3,'Color',linecolors.S); %Plot S
            semilogy(t_vals,y(:,3),'LineWidth',3,'Color',linecolors.E); %Plot E1
            semilogy(t_vals,y(:,5),'LineWidth',3,'Color',linecolors.I); %Plot I1
            semilogy(t_vals,y(:,7),'LineWidth',3,'Color',linecolors.L); %Plot L1
            semilogy(t_vals,y(:,9),'LineWidth',3,'Color',linecolors.V); %Plot V1
            ylim([1e-4 1e10]);
            xlim([0 max(t_vals)]);
            set(gca,'YScale','log','YMinorTick','off','Box','off','XTick',0:8:params.T,'XTickLabel',[],'YTick',[10.^(-4:2:10)],'YTickLabel',[cellstr(arrayfun(@(x)sprintf("$10^{%d}$",x),-4:2:10))]);
            if(deathrateindex == 2)
                set(gca,'YTickLabel',[]);
            end
            if deathrateindex == 1
                text(5,1e9,StrategyLabels{strategyindex},'FontSize',14,'FontWeight','bold');
            end
            
            set(gca,'FontSize',14,'FontWeight','bold');
           

            nexttile(length(celldeathrates)*size(PGamma{1},1)+deathrateindex);
            semilogy(t_vals,sum(y(:,3:10),2),'LineWidth',3,'Color','k','LineStyle',linestyle.(Strategies{strategyindex}));
            hold on;
            set(gca,'FontSize',14,'FontWeight','bold');
            ylim([1e3 1e10]);
            xlim([0 max(t_vals)]);
            set(gca,'YMinorTick','off','Box','off','XTick',0:8:params.T,'YTick',10.^(-2:1:10),'YTickLabel',[]);

    end

end

%% Add legends
nexttile(1);
legend({'','R','S','E','I','L','V'},'Location',[0.320178985308409,0.779965087616113,0.071621650625187,0.107379952875023]);
legend('boxoff');

nexttile(length(celldeathrates)*(size(PGamma{1},1))+1);
legend(StrategyLabels,'Location',[0.101,0.209,0.223,0.055],'Box','off');

%% Add column titles
for ii = 1:length(celldeathrates)
nexttile(ii);
title(sprintf("Cell death rate = %.2f hr$^{-1}$",celldeathrates(ii)),'Interpreter','latex','FontSize',16);
end
%% Add axis labels
xlabel(tiles,"Time (hr)",'interpreter','latex','FontSize',16);
nexttile(length(celldeathrates)+1);
ylabel('Density (mL$^{-1}$)','interpreter','latex','FontSize',16);


nexttile(length(celldeathrates)*(size(PGamma{1},1))+1);
ylabel({'Total viral genome', 'density (mL$^{-1}$)'},'interpreter','latex','FontSize',16);
set(gca,'YMinorTick','off','Box','off','XTick',0:8:params.T,'YTick',10.^(-2:1:10),'YTickLabel',arrayfun(@(x) sprintf('$10^{%d}$', x),-2:1:10,'UniformOutput',false),'TickLabelInterpreter','latex');

%% Label panels
panellabels = char('A'+(0:length(celldeathrates)*(length(PGamma))+3));

for ii = 1:length(panellabels)

    nexttile(ii);
    if mod(ii,2) == 1
        text(-6*params.T/24,10^10,sprintf("(%s)",panellabels(ii)),'FontSize',18,'FontWeight','bold');
    else
        text(-2.8*params.T/24,10^10,sprintf("(%s)",panellabels(ii)),'FontSize',18,'FontWeight','bold');
    end
end


%% Save Figure
 if SaveFigureFlag == 1

        filename = sprintf('../RevisedFigures/Figure2_TwoPanel.eps');

        exportgraphics(h,filename,'BackgroundColor','white','ContentType','vector');%save as eps
        saveas(h,[filename(1:end-4) '.svg']);%save as svg
        exportgraphics(h,[filename(1:end-4) '.png'],'BackgroundColor','white');%save as png

    end