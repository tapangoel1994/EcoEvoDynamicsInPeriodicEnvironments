%% Script generates steady state viral genome density plots as function of Q
%% for the value of gamma over which the invasion dynamics have been calculated.
%% this is to check that the excluded regions in the invasion analysis make sense. 


%% Author: Tapan Goel
%% Date: 7/18/2024


close all; 
clear all;

addpath('utils\');
addpath('lib\');
colorpalette;
fixedparameters;


CyclePeriodList = [12,16,24];
Gamma = logspace(-3,0,51);
P = linspace(0,1,51);
NumNodes = 12;
q_LV = [.2 0;.2 0;.2 0];
S0 = 1e7;
V0a = 1e4;


for index = 1:3
    if isfile(sprintf("..\\Data\\SteadyState_CyclePeriod=%.1f,S0=%1.e,V0=%1.e,p_L=%.1f,p_V=%.1f.mat",CyclePeriodList(index),S0,V0a,q_LV(index,1),q_LV(index,2)))
        load(sprintf("..\\Data\\SteadyState_CyclePeriod=%.1f,S0=%1.e,V0=%1.e,p_L=%.1f,p_V=%.1f.mat",CyclePeriodList(index),S0,V0a,q_LV(index,1),q_LV(index,2)));
        SteadyState{index} = SteadyStateDensity;
        CyclesSteadyState{index} = SSCycles;
    else
        [SteadyState{index}, CyclesSteadyState{index}] = PopulationSteadyStateFunction(CyclePeriodList(index),q_LV(index,1),q_LV(index,2),Gamma,P,NumNodes,1,params);
    end
end

% Plot the steady state viral genome densities of resident along the fixed
% gamma line

h = figure('Renderer','painters','Position',[760 200 1225*3/2 951/2]);
t = tiledlayout(1,3,'Padding','loose','TileSpacing','loose');

inset_position = [0.1150    0.335    0.18    0.4000; .425 .335 .18 .4; .713 .335 .18 .4];
colorlabel_position = [2.9109 1e4 0; 2.9109 4e3 0; 2.9109 4e2 0];
PlotTitles = {"Cycle period = 12hr","Cycle period = 16hr","Cycle period = 24hr"};
for index = 1:3
    
    SteadyState_temp = SteadyState{index};
    SteadyState_temp = squeeze(sum(SteadyState_temp(:,:,3:10),3));
    
    M = max(SteadyState_temp,[],"all","linear");
    [r,c] = find(SteadyState_temp == M);
    j = max(r);
    i = max(c);
    
    nexttile(t,index);
    semilogy(P,SteadyState_temp(:,i),'LineWidth',2,'Color','k','Marker','o','MarkerFaceColor','k','MarkerSize',4);
    hold on;
    yline(100*criticaldensitythreshold,'Color',[.5 .5 .5],'LineWidth',1.5,'LineStyle','--');
    ylim([1e-2 10*max(SteadyState_temp(:,i))]);

    set(gca,'YMinorTick','off','Box','off','XTick',0:.25:1,'YTick',logspace(-2,10,7),'XColor','r','LineWidth',2);
    set(gca,'FontSize',16,'FontWeight','bold','FontName','Times');
    title(PlotTitles{index},'FontSize',20);
    hold off;

    %% Add inset figure
    ax_inset{index} = axes('Position',inset_position(index,:));
    imagesc(P,Gamma,SteadyState_temp');    
    hold on;
    %plot(Gamma(i),Q(j),'*k','MarkerSize',5,'LineWidth',2);
    contour(P,Gamma,SteadyState_temp','k');
    yline(Gamma(i),'LineWidth',1.5,'Color','r');
    yticks([1e-3 1e-2 1e-1 1e0]);
    yticklabels({'10$^{-3}$','10$^{-2}$','10$^{-1}$','1'});
    xticks(linspace(0,1,5));
    c = colorbar;
    c.Label.String = "Viral genome density (mL$^{-1}$)";
    c.Label.Interpreter = 'latex';
    c.Label.Rotation = -90;
    c.Label.Position = colorlabel_position(index,:);
    
    set(ax_inset{index},'YScale','log','XScale','linear','ColorScale','log','YDir','normal'); 
    ylabel('$\gamma$: Induction rate (hr$^{-1}$)');
    xlabel( '$p$: Integration Probability');
    pbaspect(ax_inset{index},[1 1 1]);
    hold off;

end

xlabel(t,'$p$: Integration probability','FontSize',18,'Interpreter','latex');
ylabel(t,{'Viral genome density (mL$^{-1}$)'},'FontSize',18,'Interpreter','latex');

%% Add plot labels
nexttile(t,1);
text(-.21,1e8,'(A)','FontSize',20,'FontWeight','bold');
text(1.13,1e8,'(B)','FontSize',20,'FontWeight','bold');
text(2.42,1e8,'(C)','FontSize',20,'FontWeight','bold');

%% Save Figure
filename = dir('..\\Figures\\FigureS6*');

if isempty(filename)
    filename = '..\\Figures\\FigureS6_v1.eps';
else
    filename = [filename(end).folder '\' filename(end).name];
    version = extractBetween(filename,"_v",".");
    version = version{1};
    version = str2num(version);
    filename = [extractBefore(filename,['v' num2str(version)]) 'v' num2str(version+1) '.eps'];
end
filename1 = [filename(1:end-4) '.png'];

exportgraphics(h,filename,"BackgroundColor",'none','ContentType','vector');
exportgraphics(h,filename1,"BackgroundColor",'white');


   

