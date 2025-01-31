%% Code generates figure 3 and figure S5
%% Date created: 6/13/2024
%% Author: Tapan Goel

close all; 
clear all;

%% Figure making conditions
SaveFigureFlag = 0; % if SaveFigureFlag is 1, save the figure as eps and as png.
PlotAllDeathRates = 1; %if this variable is 1, generate figures for both low and high death rates, %if zero, generate figure for only the low death rate case
%% Add file paths
addpath('utils/','lib/','../Data/');
colorpalette;
VariableMaps;

%% Set Cycle period and filtration conditions
Q_LV = [0 .2;.1 .1;.2 0]; %each row contains q_L and q_V for a particular filtration condition
CyclePeriod = 24; % growth cycle period in hours

%% Set strategies and labels

StrategyLabels = {'Obligately lytic','Temperate','Obligately lysogenic'};
ParameterFiles = {'lowdeathrate','highdeathrate'};
Filtrate = {'virions only','virions and lysogens','lysogens only'};

%% Create maps:
 filtertotile = [1 3;6 8;11 13]; %filtertotile(filtrationindex,1) maps to the tiles corresponding to the plot for that particular filtration condition
 plotlabel = ['A' 'B';'C' 'D';'E' 'F'];



for CellDeathRatesIndex = 1:(PlotAllDeathRates+1) % separate figure for each death rate
        
        %% Load parametervalues
        run(['fixedparameters_' ParameterFiles{CellDeathRatesIndex}]);
        %% Simulation parameters
        params.T = CyclePeriod;
        params.t_vals = transpose(0:params.dt:params.T);
        %% initial conditions
        R0 = 100; %initial resource amount in ug/mL ( 500 mL flask)
        S0 = 1e7; %Initial concentration of susceptibles in flask (per mL)
        Va_0= 1e4; %initial concentration of virus in flask (per mL)
        Vb_0 = 0;
        x0 = [R0 S0 zeros(1,6) Va_0 Vb_0];


        h_main = figure('Position',[100 0 1400 800]);
        t_main = tiledlayout(3,5,'TileSpacing','loose','Padding','loose');
        
              
        NCycles = 40;

        for filtrationindex = 1:size(Q_LV,1) %% separate row for each filtration condition
            q_L = Q_LV(filtrationindex,1);
            q_V = Q_LV(filtrationindex,2);
            
            TransferMatrix = diag([q_R q_S q_E q_E q_I q_I q_L q_L q_V q_V]);

            for strategyindex = 1:length(Strategies)

                params.p = [PGamma.(Strategies{strategyindex})(1) 0];
                params.gamma = [PGamma.(Strategies{strategyindex})(2) 0];
                
                T = [];
                Y = [];
                x0 = [R0 S0 zeros(1,6) Va_0 Vb_0];
                
                for iter = 1:NCycles   %% simulate 40 cycles
                    [t_vals, y] = ode113(@ODE_RSEILV_2Species, params.t_vals, x0, options, params);
                    T = [T; params.T*(iter-1) + t_vals];
                    Y = [Y;y];
                    x0 = [R0 S0 zeros(1,8)] + y(end,:)*TransferMatrix;
                end
                
                ViralGenomeDensity = sum(Y(:,3:10),2); % calculate total density of viral genomes
                ViralGenomeDensity(ViralGenomeDensity < criticaldensitythreshold/10) = criticaldensitythreshold/10; %artificially set viral genome density to 1/10 of the critical threshold to be able to plot lines around it
                CycleStartDensities = [T(1:length(t_vals):end) ViralGenomeDensity(1:length(t_vals):end)]; % store total viral genome density at beginning of a cycle

                %definitely plot the viral density for the first 4 cycles
                %and after that, set the viral density to zero if its
                %continuous zeros
                x = find(CycleStartDensities(:,2) == criticaldensitythreshold/10);
                x(x <= 7) = [];
                CycleStartDensities(x,2) = 0;

                figure(h_main);

                %% Generate left panel with the dynamics of the first 3 cycles
                nexttile(t_main,filtertotile(filtrationindex,1),[1 2]);

                semilogy(T(1:3*length(t_vals)+1),ViralGenomeDensity(1:3*length(t_vals)+1),'LineWidth',1.5,'Color','k','LineStyle',getfield(linestyle,Strategies{strategyindex}));
                hold on;
                semilogy(T(1:length(t_vals):3*length(t_vals)+1),ViralGenomeDensity(1:length(t_vals):3*length(t_vals)+1),'Marker',getfield(marker,Strategies{strategyindex}),'Color','k','LineStyle','none');
                
                %% Generate right panel for cycle-to-cycle dynamics
                nexttile(t_main,filtertotile(filtrationindex,2),[1 3]);

                semilogy(CycleStartDensities(:,1),CycleStartDensities(:,2),'LineWidth',1,'Color','k','LineStyle',getfield(linestyle,Strategies{strategyindex}));
                hold on;
                semilogy(CycleStartDensities([1:4 6:2:end],1),CycleStartDensities([1:4 6:2:end],2),'Marker',getfield(marker,Strategies{strategyindex}),'Color','k','LineStyle','none');
                
            
            end
            
            %% Annotate and format left plots
            nexttile(t_main,filtertotile(filtrationindex,1),[1 2]);
            
            xline(0.99*params.T,'--','LineWidth',1.5,'Color',.8*[1 1 1]);
            xline(1.99*params.T,'--','LineWidth',1.5,'Color',.5*[1 1 1]);
            xline(2.99*params.T,'--','LineWidth',1.5,'Color',.2*[1 1 1]);
            yline(criticaldensitythreshold,':r','LineWidth',1);
            %area([0 3.1*params.T],[1e-2 1e-2],'FaceColor',[.5 .5 .5],'FaceAlpha',1);
            ylim([1e-4 1e10]);
            xlim([0 3.1*params.T]);
            set(gca,'YMinorTick','off','Box','off','XTick',0:params.T/3:3*params.T,'XTickLabel',[],'YTick',10.^(-4:2:10),'YTickLabel',arrayfun(@(x)sprintf("$10^{%d}$",x),-4:2:10));
            set(gca,'FontSize',14,'FontWeight','bold','FontName','Times');
            text(-12*params.T/24,2e10,sprintf('(%s)',plotlabel(filtrationindex,1)),'FontSize',14,'FontWeight','bold');
            text(.3*params.T,1e-1,'Cycle 1','FontSize',14,'FontWeight','bold');
            text(1.3*params.T,1e-1,'Cycle 2','FontSize',14,'FontWeight','bold');
            text(2.3*params.T,1e-1,'Cycle 3','FontSize',14,'FontWeight','bold');
            title(sprintf("Filtrate: %s",Filtrate{filtrationindex}),'FontSize',16,'FontWeight','bold','FontName','Times','Interpreter','latex');
            hold off;
            
            %% Annotate and format right plot
            nexttile(t_main,filtertotile(filtrationindex,2),[1 3]);
            
            xline(0.99*params.T,'--','LineWidth',1.5,'Color',.8*[1 1 1]);
            xline(1.99*params.T,'--','LineWidth',1.5,'Color',.5*[1 1 1]);
            xline(2.99*params.T,'--','LineWidth',1.5,'Color',.2*[1 1 1]);
            yline(criticaldensitythreshold,':r','LineWidth',1);
            %area([0 (NCycles)*params.T],[1e-2 1e-2],'FaceColor',[.5 .5 .5],'FaceAlpha',.95);
            ylim([1e-4 1e10]);
            xlim([0 (NCycles)*params.T]);
            set(gca,'YMinorTick','off','Box','off','XTick',[0:params.T:3*params.T 5*params.T:2*params.T:T(end)],'XTickLabel',[],'YTick',10.^(-4:2:10),'YTickLabel',[]);
            set(gca,'FontSize',14,'FontWeight','bold','FontName','Times');
            text(-50*params.T/24,2e10,sprintf('(%s)',plotlabel(filtrationindex,2)),'FontSize',16,'FontWeight','bold');
            title(sprintf("Filtrate: %s",Filtrate{filtrationindex}),'FontSize',14,'FontWeight','bold','FontName','Times','Interpreter','latex');
            hold off;
            

        end
        
        %% Set tick labels, axes labels and legends
        nexttile(t_main,11,[1 2]);
        set(gca,'XTick',0:params.T/3:3*params.T,'XTickLabel',0:params.T/3:3*params.T);
        
        xlabel('Time (hr)','FontSize',16,'FontWeight','bold','FontName','Times');

        nexttile(t_main,13,[1 3]);
        set(gca,'XTick',[0:params.T:3*params.T 5*params.T:2*params.T:T(end)],'XTickLabel',[0:1:3 5:2:NCycles]);
        xlabel('Number of cycles','FontSize',16,'FontWeight','bold','FontName','Times');
        
        nexttile(t_main,6,[1 2]);
        ylabel('Total viral genome density (mL$^{-1}$)','FontSize',16,'FontWeight','bold','FontName','Times','Interpreter','latex','Position',[-12,2236.1,-1]);

        nexttile(t_main,3,[1 3]);
        hold on;
        lh(1) = plot(nan,nan,'LineStyle',linestyle.lytic,'Marker',marker.lytic,'LineWidth',1,'Color','k');
        lh(2) = plot(nan,nan,'LineStyle',linestyle.temperate,'Marker',marker.temperate,'LineWidth',1,'Color','k');
        lh(3) = plot(nan,nan,'LineStyle',linestyle.lysogenic,'Marker',marker.lysogenic,'LineWidth',1,'Color','k');
        legend(lh,Strategies,'Box','off','Location',[0.79,0.74,0.089,0.082]);


        %% Save figure
        if SaveFigureFlag == 1
            filename = sprintf('../RevisedFigures/Figure3_%s.eps',ParameterFiles{CellDeathRatesIndex});
            saveas(h_main,[filename(1:end-4) '.svg']);%save as svg
        end

end
        