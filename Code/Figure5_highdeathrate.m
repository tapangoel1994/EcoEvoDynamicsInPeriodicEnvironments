%% This script generates Figure S11 for the revised manuscript.
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
fixedparameters_highdeathrate;

%% Define parameters
Gamma = 10.^(-6:.1:0);
P = 0:.02:1;


CyclePeriodList = [4, 8, 12, 16, 20, 24, 36, 48, 96];

q_LV = [0.025 0;0.05 0;0.1 0;0.2 0;0.4 0; ...
    0 .05; 0 .1; 0 .2];

filtratetype = {'lysogens only','virions only'};

columnnames = {'CyclePeriod','q_L','q_V','PMax','GammaMax','MaxDensity','SteadyStateHeatMap','SteadyStateCycles','pESS','InvasionSuccessMatrix'};
    
OptimalStrategy = table('Size',[0,length(columnnames)],'VariableTypes',[repmat("double",1,6),"cell","cell","double","cell"],'VariableNames',columnnames);

SteadyState = cell(length(CyclePeriodList),size(q_LV,1));
CyclesSteadyState = cell(length(CyclePeriodList),size(q_LV,1));
Invasion = cell(length(CyclePeriodList),size(q_LV,1));
CyclesInvasion = cell(length(CyclePeriodList),size(q_LV,1));
InvasionSuccessMatrix = cell(length(CyclePeriodList),size(q_LV,1));

for cycleperiodindex = 1:length(CyclePeriodList)
    for filterindex = 1:size(q_LV,1)

        %% Load steady state density matrices

        filename = sprintf('../Data/SteadyState_CyclePeriod=%d,d=%.2f,q_L=%.3f,q_V=%.3f.mat',CyclePeriodList(cycleperiodindex), ...
            params.d_S,q_LV(filterindex,1),q_LV(filterindex,2));
        if isfile(filename)
            load(filename);
            SteadyState{cycleperiodindex,filterindex} = SteadyStateDensity;
            CyclesSteadyState{cycleperiodindex,filterindex} = SSCycles;
            
        else
            disp(sprintf("Steady State file does not exist for T = %d, q_L = %.3f, q_V = %.3f",CyclePeriodList(cycleperiodindex),q_LV(filterindex,1),q_LV(filterindex,2)));
        end
        
        %% Find population dynamics optimal strategy
        SteadyState_temp = SteadyState{cycleperiodindex,filterindex};
        SteadyState_temp = squeeze(sum(SteadyState_temp(:,:,3:10),3));
        
        [r,c] = find(SteadyState_temp == max(SteadyState_temp,[],"all","linear"));
        j = max(r);
        i = max(c);

        if max(SteadyState_temp,[],"all","linear") > 10*criticaldensitythreshold % doing invasion analysis for the optimal strategy only makes sense if the strategy produces more viral genomes than the lower threshold
            
            
            filename = sprintf('../Data/Invasion_CyclePeriod=%d,d=%.2f,Gamma=%.3e,q_L=%.3f,q_V=%.3f.mat',CyclePeriodList(cycleperiodindex), ...
                params.d_S,Gamma(i),q_LV(filterindex,1),q_LV(filterindex,2));

            if isfile(filename)
                
                load(filename,"InvasionDensity", "InvasionMatrix", "CyclesToInvasion");
                Invasion{cycleperiodindex,filterindex} = InvasionDensity;
                CyclesInvasion{cycleperiodindex,filterindex} = CyclesToInvasion;
                InvasionSuccessMatrix{cycleperiodindex,filterindex} = InvasionMatrix;
                
            else
                disp(sprintf("Invasion file does not exist for T = %d, q_L = %.3f, q_V = %.3f",CyclePeriodList(cycleperiodindex),q_LV(filterindex,1),q_LV(filterindex,2)));
            end

       
        
        %% Find ESS Strategy
                if ~isempty(InvasionSuccessMatrix{cycleperiodindex,filterindex})
                    x = FindESSStrategy(InvasionSuccessMatrix{cycleperiodindex,filterindex}',P);
                    if ~isempty(x)
                        
                        
                        temp = {CyclePeriodList(cycleperiodindex), q_LV(filterindex,1), q_LV(filterindex,2),P(j),Gamma(i),max(SteadyState_temp,[],"all","linear") ...
                            , SteadyState(cycleperiodindex,filterindex), CyclesSteadyState(cycleperiodindex,filterindex), x, InvasionSuccessMatrix(cycleperiodindex,filterindex)};
                        
                        OptimalStrategy = [OptimalStrategy;cell2table(temp,VariableNames=columnnames)];
                    else
                        temp = {CyclePeriodList(cycleperiodindex), q_LV(filterindex,1), q_LV(filterindex,2),P(j),Gamma(i),max(SteadyState_temp,[],"all","linear") ...
                            , SteadyState(cycleperiodindex,filterindex), CyclesSteadyState(cycleperiodindex,filterindex), NaN(1), InvasionSuccessMatrix(cycleperiodindex,filterindex)};
                        
                        OptimalStrategy = [OptimalStrategy;cell2table(temp,VariableNames=columnnames)];

                    end
                end
        end

    end
end


%% Generate figure;

h = figure('Renderer','painters','Position',[318,55,1196,1275]);
t = tiledlayout(3,2,"TileSpacing","loose","Padding","loose");

%% Plot PIPs
optimalstrategymarkers = {'^','d','s','o'};
markeridx = 1;
for qV = [.2 0]
    for cycleperiod = [16 24]

        tile = nexttile;

        currData = OptimalStrategy(OptimalStrategy.q_L == 0.2-qV & OptimalStrategy.q_V == qV & OptimalStrategy.CyclePeriod == cycleperiod,:);

        imagesc(P,P,currData.InvasionSuccessMatrix{:}');
        colormap(tile,PIPColorMap); 
        hold on;
        plot(currData.pESS,currData.pESS,'MarkerSize',5,'LineWidth',2,'LineStyle','none','Marker',optimalstrategymarkers{markeridx},'MarkerFaceColor','r','MarkerEdgeColor','r'); markeridx = markeridx+1;
        xticks(linspace(0,1,5));
        xticklabels(linspace(0,1,5));
        yticks(linspace(0,1,5));
        yticklabels(linspace(0,1,5));
        set(gca,'PlotBoxAspectRatio',[1 1 1]);
        set(gca,'YDir','normal','FontSize',14,'FontWeight','bold','CLim',[-2 1]);
        xlabel('p (resident)','FontSize',16,'FontWeight','bold'); 
        ylabel('p (mutant)','FontSize',16,'FontWeight','bold');
        title({sprintf('Filtrate: $q_V = %.2f,q_L = %.2f$,',qV,.2-qV), sprintf('Cycle period = %d hr',cycleperiod)},'FontSize',16,'Interpreter','latex','FontWeight','bold');

    end
end


%% plot ESS Strategies as function of time period for lysogen filtration cases

q_LV_sub = unique(OptimalStrategy(OptimalStrategy.q_V == 0,:).q_L);
q_LV_sub = [q_LV_sub zeros(size(q_LV_sub))];

for idx = 1:size(q_LV_sub,1)

    currData = OptimalStrategy(OptimalStrategy.q_V == q_LV_sub(idx,2) & OptimalStrategy.q_L == q_LV_sub(idx,1),:);
    currData = currData(~isnan(currData.pESS),:);

    nexttile(5);
    
    plot(currData.CyclePeriod,currData.pESS,'Color',idx*[1 0 0]/(size(q_LV_sub,1)+1),'Marker','o','LineWidth',1.5); %% plot the pESS
    hold on;

   
    nexttile(6);
    plot(currData.CyclePeriod, currData.PMax - currData.pESS,'Color', idx*[1 0 0]/(size(q_LV_sub,1)+1),'Marker','o','LineWidth',1.5); %% plot the gamma*s where the p* is not extreme.
    hold on;
end


nexttile(5);
    
    plot(OptimalStrategy(OptimalStrategy.q_L == .2 & OptimalStrategy.q_V == 0 & OptimalStrategy.CyclePeriod == 16,:).CyclePeriod, ...
        OptimalStrategy(OptimalStrategy.q_L == .2 & OptimalStrategy.q_V == 0 & OptimalStrategy.CyclePeriod == 16,:).pESS, ...
        'LineStyle','none','Marker','s','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',10);
    
    plot(OptimalStrategy(OptimalStrategy.q_L == .2 & OptimalStrategy.q_V == 0 & OptimalStrategy.CyclePeriod == 24,:).CyclePeriod, ...
        OptimalStrategy(OptimalStrategy.q_L == .2 & OptimalStrategy.q_V == 0 & OptimalStrategy.CyclePeriod == 24,:).pESS, ...
        'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',10);
    
    lgd = legend(arrayfun(@(x)sprintf('$q_L = %.3f$',x),q_LV_sub(:,1),'UniformOutput',false),'Interpreter','latex','Location','northeast','FontSize',14,'Box','off');
    set(gca,'XTick',0:4:24,'XTickLabel',0:4:24,'YLim',[0 1],'YTick',0:.1:1,'FontSize',14,'FontWeight','bold','TickLabelInterpreter','latex','Box','off');
    xlabel('Cycle period (hr)','FontSize',16,'Interpreter','latex','FontWeight','bold');
    ylabel('$p_{ESS}$','FontSize',16,'Interpreter','latex','FontWeight','bold','Rotation',90);
    
    
    
    nexttile(6);
    
    plot(OptimalStrategy(OptimalStrategy.q_L == .2 & OptimalStrategy.q_V == 0 & OptimalStrategy.CyclePeriod == 16,:).CyclePeriod, ...
        OptimalStrategy(OptimalStrategy.q_L == .2 & OptimalStrategy.q_V == 0 & OptimalStrategy.CyclePeriod == 16,:).PMax - ...
        OptimalStrategy(OptimalStrategy.q_L == .2 & OptimalStrategy.q_V == 0 & OptimalStrategy.CyclePeriod == 16,:).pESS, ...
        'LineStyle','none','Marker','s','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',10);
    
    plot(OptimalStrategy(OptimalStrategy.q_L == .2 & OptimalStrategy.q_V == 0 & OptimalStrategy.CyclePeriod == 24,:).CyclePeriod, ...
        OptimalStrategy(OptimalStrategy.q_L == .2 & OptimalStrategy.q_V == 0 & OptimalStrategy.CyclePeriod == 24,:).PMax - ...
        OptimalStrategy(OptimalStrategy.q_L == .2 & OptimalStrategy.q_V == 0 & OptimalStrategy.CyclePeriod == 24,:).pESS, ...
        'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',10);
    
    lgd = legend(arrayfun(@(x)sprintf('$q_L = %.3f$',x),q_LV_sub(:,1),'UniformOutput',false),'Interpreter','latex','Location','southeast','FontSize',14,'Box','off');
    set(gca,'XTick',0:4:24,'XTickLabel',0:4:24,'YLim',[0 1],'YTick',0:.1:1,'FontSize',14,'FontWeight','bold','TickLabelInterpreter','latex','Box','off', 'YMinorTick','off');
    xlabel('Cycle period (hr)','FontSize',16,'Interpreter','latex','FontWeight','bold');
    ylabel('$p^* - p_{ESS}$','FontSize',16,'Interpreter','latex','FontWeight','bold','Rotation',90);
    

panellabels = 'ABCDEF';

for index = 1:length(panellabels)
    
    nexttile(index);
    
    xloc = get(gca,'XTick');

    xloc = min(xloc);

    yloc = get(gca,'YTick');

    yloc = max(yloc);

    text(xloc,yloc, sprintf("(%s)",panellabels(index)),'FontSize',18,'FontWeight','bold');

end

%% save figure

saveas(h,'../RevisedFigures/FigureS11.svg');