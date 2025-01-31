%% This script generates Figure S10 for the revised manuscript.
%For a range of cycle periods and filtration conditions, this script
%generates the steady state density plots as functions of viral strategies.

%% Note: Run script DataGenerate_highdeathrate before running this script.
%% Note: Figure has annotations on it that need to be moved around in inkscape for better visibility. 

%% Author: Tapan Goel
%% Date last modified: 1/31/2025

clear all;
close all;

addpath('utils/','lib/','../Data/');
colorpalette;
VariableMaps;

for deathrateidx = 2 %% run script for high deathrate case

    run(['fixedparameters_' FixedParameterSets{deathrateidx}]);

    %% Define parameters
    if strcmp('highdeathrate',FixedParameterSets{deathrateidx})
        Gamma = 10.^(-6:.1:0);
    else
        Gamma = 10.^(-6:.1:-1);
    end

    P = 0:.02:1;
    

    NumNodes = 64;
    
    
    CyclePeriodList = [4, 8, 12, 16, 20, 24, 36, 48, 96];
    
    q_LV = [0.025 0;0.05 0;0.1 0;0.2 0;0.4 0;...
            0 0.05;0 0.1;0 0.2];
    
    columnnames = {'CyclePeriod','q_L','q_V','PMax','GammaMax','MaxDensity','SteadyStateHeatMap','SteadyStateCycles'};
    
    OptimalStrategy = table('Size',[0,length(columnnames)],'VariableTypes',[repmat("double",1,length(columnnames)-2),"cell","cell"],'VariableNames',columnnames);
    
    %% Load steady state density matrices from files. If file doesnt exist, give an error message.
    for cycleperiodindex = 1:length(CyclePeriodList)
        for filterindex = 1:size(q_LV,1)
            filename = sprintf('../Data/SteadyState_CyclePeriod=%d,d=%.2f,q_L=%.3f,q_V=%.3f.mat',CyclePeriodList(cycleperiodindex), ...
                params.d_S,q_LV(filterindex,1),q_LV(filterindex,2));
            if isfile(filename)
                load(filename,"SteadyStateDensity","SSCycles");
                SteadyState{cycleperiodindex,filterindex} = SteadyStateDensity;
                CyclesSteadyState{cycleperiodindex,filterindex} = SSCycles;
    
                SteadyState_temp = SteadyState{cycleperiodindex,filterindex};
                SteadyState_temp = squeeze(sum(SteadyState_temp(:,:,3:10),3));
                
                [r,c] = find(SteadyState_temp == max(SteadyState_temp,[],"all","linear"));
                j = max(r);
                i = max(c);
                
                newrow = {CyclePeriodList(cycleperiodindex), q_LV(filterindex,1), q_LV(filterindex,2), P(j) Gamma(i), max(SteadyState_temp,[],"all","linear"),{SteadyState_temp},{SSCycles}};
                OptimalStrategy = [OptimalStrategy; ...
                cell2table(newrow,'VariableNames',columnnames) ];
            else
                disp(sprintf("File does not exist for T = %d, q_L = %.3f, q_V = %.3f",CyclePeriodList(cycleperiodindex),q_LV(filterindex,1),q_LV(filterindex,2)));
                Save.Flag = 1;
                Save.FileName = filename;
                [SteadyState{cycleperiodindex,filterindex}, CyclesSteadyState{cycleperiodindex,filterindex}] = PopulationSteadyStateFunction(CyclePeriodList(cycleperiodindex), ...
                   q_LV(filterindex,1),q_LV(filterindex,2),MaxCycles,Gamma,P,12,Save,params);
            end
    
            
        end
    end
    
    
    h = figure('Renderer','painters','Position',[431,82,979,1210]);
    t = tiledlayout(3,2,'Padding','loose','TileSpacing','loose');
    
    
    %% Plot HeatMaps.
    optimalstrategymarkers = {'^','d','s','o'};
    markeridx = 1;
    for qV = [.2 0]
        for cycleperiod = [16 24]
        
    
            nexttile;
            
            currData = OptimalStrategy(OptimalStrategy.q_L == .2-qV & OptimalStrategy.q_V == qV & OptimalStrategy.CyclePeriod == cycleperiod,:);
            TotalSteadyStateDensity = currData.SteadyStateHeatMap{1};
    
            imagesc(P,Gamma,TotalSteadyStateDensity');
            hold on;
            plot(currData.PMax,currData.GammaMax,'MarkerSize',5,'LineWidth',2,'LineStyle','none','Marker',optimalstrategymarkers{markeridx},'MarkerFaceColor','r','MarkerEdgeColor','r'); markeridx = markeridx+1;
            
            contour(P,Gamma,TotalSteadyStateDensity','k');
            yticks(10.^(-6:1:0));
            yticklabels(arrayfun(@(x)sprintf("10$^{%d}$",x),-6:1:0));
            xticks(linspace(0,1,5));
            c = colorbar;
            c.Label.String = "Viral genome density (mL$^{-1}$)";
            c.Label.Interpreter = 'latex';
            c.Label.Rotation = -90;
            c.FontSize = 14;
            c.Label.FontSize = 16;
            
            set(gca,'YScale','log','XScale','linear','ColorScale','log','YDir','normal','FontWeight','bold','FontSize',14); 
            ylabel('$\gamma$: Induction rate (hr$^{-1}$)','FontSize',16,"FontWeight",'bold','Interpreter','latex');
            xlabel('$p$: Integration probability','FontSize',16,"FontWeight",'bold','Interpreter','latex');
            
            text(.60,.0295,'$(p^*,\gamma^*)$','interpreter','latex','FontSize',14,'FontWeight','bold','HorizontalAlignment','right','Color','k','BackgroundColor','white');
            text(.52,.00832,'$\longrightarrow$','interpreter','latex','FontSize',14,'FontWeight','bold','HorizontalAlignment','right','Color','k','Rotation',-90);
            
            title({sprintf('Filtrate: $q_V = %.2f,q_L = %.2f$,',qV,.2-qV), sprintf('Cycle period = %d hr',cycleperiod)},'FontSize',16,'Interpreter','latex','FontWeight','bold');
            c.Label.Position(1) = 4.7;
        end
    end
    
    
    %% Plot optimal strategies for only lysogen passage
    
    q_LV_sub = q_LV(q_LV(:,2)==0,:);
    
    for idx = 1:size(q_LV_sub,1)
    
        currData = OptimalStrategy(OptimalStrategy.q_L == q_LV_sub(idx,1) & OptimalStrategy.q_V == q_LV_sub(idx,2) & OptimalStrategy.MaxDensity > criticaldensitythreshold,:);
        
        nexttile(5);
    
        plot(currData.CyclePeriod,currData.PMax,'Color',idx*[1 1 1]/(size(q_LV_sub,1)+1),'Marker','o','LineWidth',1.5); %% plot the p*'s
        hold on;
    
        currData_sub = currData(currData.GammaMax > min(Gamma) & currData.GammaMax < max(Gamma),:);
        
        nexttile(6);
        semilogy(currData_sub.CyclePeriod, currData_sub.GammaMax,'Color', idx*[1 1 1]/(size(q_LV_sub,1)+1),'Marker','o','LineWidth',1.5); %% plot the gamma*s where the p* is not extreme.
        hold on;
    
    end
    
    nexttile(5);
    
    plot(OptimalStrategy(OptimalStrategy.q_L == .2 & OptimalStrategy.q_V == 0 & OptimalStrategy.CyclePeriod == 24,:).CyclePeriod, ...
        OptimalStrategy(OptimalStrategy.q_L == .2 & OptimalStrategy.q_V == 0 & OptimalStrategy.CyclePeriod == 24,:).PMax, ...
        'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',10);
    
    plot(OptimalStrategy(OptimalStrategy.q_L == .2 & OptimalStrategy.q_V == 0 & OptimalStrategy.CyclePeriod == 16,:).CyclePeriod, ...
        OptimalStrategy(OptimalStrategy.q_L == .2 & OptimalStrategy.q_V == 0 & OptimalStrategy.CyclePeriod == 16,:).PMax, ...
        'LineStyle','none','Marker','s','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',10);
    
    lgd = legend(arrayfun(@(x)sprintf('$q_L = %.3f$',x),q_LV_sub(:,1),'UniformOutput',false),'Interpreter','latex','Location','southeast','FontSize',14,'Box','off');
    set(gca,'XTick',0:4:24,'XTickLabel',0:4:24,'YLim',[0 1],'YTick',0:.1:1,'FontSize',14,'FontWeight','bold','TickLabelInterpreter','latex','Box','off');
    xlabel('Cycle period (hr)','FontSize',16,'Interpreter','latex','FontWeight','bold');
    ylabel('$p^*$','FontSize',16,'Interpreter','latex','FontWeight','bold','Rotation',0);
    
    
    
    nexttile(6);
    
    semilogy(OptimalStrategy(OptimalStrategy.q_L == .2 & OptimalStrategy.q_V == 0 & OptimalStrategy.CyclePeriod == 24,:).CyclePeriod, ...
        OptimalStrategy(OptimalStrategy.q_L == .2 & OptimalStrategy.q_V == 0 & OptimalStrategy.CyclePeriod == 24,:).GammaMax, ...
        'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',10);
    
    semilogy(OptimalStrategy(OptimalStrategy.q_L == .2 & OptimalStrategy.q_V == 0 & OptimalStrategy.CyclePeriod == 16,:).CyclePeriod, ...
        OptimalStrategy(OptimalStrategy.q_L == .2 & OptimalStrategy.q_V == 0 & OptimalStrategy.CyclePeriod == 16,:).GammaMax, ...
        'LineStyle','none','Marker','s','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',10);
    
    lgd = legend(arrayfun(@(x)sprintf('$q_L = %.3f$',x),q_LV_sub(:,1),'UniformOutput',false),'Interpreter','latex','Location','southeast','FontSize',14,'Box','off');
    set(gca,'XLim',[0 24],'XTick',0:4:24,'XTickLabel',0:4:24,'YLim',[5e-7 1e0],'YTick',10.^(-6:1:0),'FontSize',14,'FontWeight','bold','TickLabelInterpreter','latex','Box','off', 'YMinorTick','off');
    xlabel('Cycle period (hr)','FontSize',16,'Interpreter','latex','FontWeight','bold');
    ylabel('$\gamma^*$','FontSize',16,'Interpreter','latex','FontWeight','bold','Rotation',0);
    
    
    
    %% Add Panel Labels
    
    panellabels = 'ABCDEF';
    
    for index = 1:length(panellabels)
    
        nexttile(index);
    
        xloc = get(gca,'XTick');
    
        xloc = min(xloc)- .4*(max(xloc)-min(xloc));
    
        yloc = get(gca,'YTick');
    
        yloc = max(yloc);
    
        text(xloc,yloc, sprintf("(%s)",panellabels(index)),'FontSize',18,'FontWeight','bold');
    end

end



%% Save figure
saveas(h,'../RevisedFigures/FigureS10.svg');