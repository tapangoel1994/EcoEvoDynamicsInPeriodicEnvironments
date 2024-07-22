%%Code to generate figure 5 of the manuscript: Steady State Density and PIP
%%plots for passaging only virions and only lysogens.

%% Date created: 05/29/2024
%% Author: Tapan Goel

%%Code pulled from scripts for steady state analysis and PIP analysis using
%%parallelization

close all; 
clear all;

addpath('utils\');
addpath('lib\');
colorpalette;
fixedparameters;

CyclePeriodList = [12, 16,24];
Gamma = logspace(-3,0,51);
Q = linspace(0,1,51);
NumNodes = 12;
p_LV = [.2 0;.2 0;.2 0];
S0 = 1e7;
V01 = 1e4;


for index = 1:3
    if isfile(sprintf("..\\Data\\SteadyState_CyclePeriod=%.1f,S0=%1.e,V0=%1.e,p_L=%.1f,p_V=%.1f.mat",CyclePeriodList(index),S0,V01,p_LV(1,1),p_LV(1,2)))
        load(sprintf("..\\Data\\SteadyState_CyclePeriod=%.1f,S0=%1.e,V0=%1.e,p_L=%.1f,p_V=%.1f.mat",CyclePeriodList(index),S0,V01,p_LV(1,1),p_LV(1,2)));
        SteadyState{index} = SteadyStateDensity;
        CyclesSteadyState{index} = SSCycles;
    else
        [SteadyState{index}, CyclesSteadyState{index}] = PopulationSteadyStateFunction(CyclePeriodList(index),p_LV(1,1),p_LV(1,2),Gamma,Q,NumNodes,1,params);
    end
end

for index = 1:3
    
    SteadyState_temp = SteadyState{index};
    [M,I] = max(SteadyState_temp(:,:,7),[],"all","linear");
    [j,i]=ind2sub([length(Q),length(Gamma)],I);
    
    
    InvasionVariable = [Q' Gamma(i)*ones(size(Q'))];
   
    if isfile(sprintf("..\\Data\\Invasion_CyclePeriod=%.1f,S0=%1.e,V0=%1.e,p_L=%.1f,p_V=%.1f.mat",CyclePeriodList(index),S0,Va_0,p_LV(index,1),p_LV(index,2)))
        load(sprintf("..\\Data\\Invasion_CyclePeriod=%.1f,S0=%1.e,V0=%1.e,p_L=%.1f,p_V=%.1f.mat",CyclePeriodList(index),S0,Va_0,p_LV(index,1),p_LV(index,2)));
        Invasion{index} = InvasionDensity;
        CyclesInvasion{index} = CyclesToInvasion;
        InvasionSuccessMatrix{index} = InvasionMatrix;

    else
        [Invasion{index}, InvasionSuccessMatrix{index}, CyclesInvasion{index}] = InvasionDynamics(CyclePeriodList(index),p_LV(index,1),p_LV(index,2),InvasionVariable,NumNodes,1,params);
        
    end


end

MaximaPoint = {};
PlotTitles = {'Cycle period = 12hr','Cycle period = 16hr','Cycle period = 24hr'};

%% Plot 
h = figure('Renderer','painters','Position',[760 10 1225 .9*1426]);
t = tiledlayout(3,2,'Padding','loose','TileSpacing','loose');

for index = 1:3
    
    SteadyState_temp = SteadyState{index};
    SteadyState_temp = squeeze(sum(SteadyState_temp(:,:,3:10),3));
    InvasionSuccess_temp = InvasionSuccessMatrix{index};

    M = max(SteadyState_temp,[],"all","linear");
    [r,c] = find(SteadyState_temp == M);
    j = max(r);
    i = max(c);
    MaximaPoint{index} = [Gamma(i) Q(j)];

    
    
    %% Plot steady state densities
    tile = nexttile();
    imagesc(Q,Gamma,SteadyState_temp');
    hold on;
    plot(Q(j),Gamma(i),'*k','MarkerSize',5,'LineWidth',2);
    
    
    
    contour(Q,Gamma,SteadyState_temp','k');
    yticks([1e-3 1e-2 1e-1 1e0]);
    yticklabels({'10$^{-3}$','10$^{-2}$','10$^{-1}$','1'});
    xticks(linspace(0,1,5));
    c = colorbar;
    c.Label.String = "Viral genome density (mL$^{-1}$)";
    c.Label.Interpreter = 'latex';
    c.Label.Rotation = -90;
    c.FontSize = 12;
    c.Label.FontSize = 14;
    
    %yticklabels({'0','.5','1'});
    set(gca,'YScale','log','XScale','linear','ColorScale','log','YDir','normal','FontWeight','bold','FontSize',12); 
    ylabel('$\gamma$: Induction rate (hr$^{-1}$)','FontSize',14,"FontWeight",'bold');
    xlabel( '$p$: Integration Probability','FontSize',14,"FontWeight",'bold');
    c.Label.Position(1) = 4.3167;
    hold off;

    %% Plot PIP
    tile = nexttile;
    imagesc(Q,Q,InvasionSuccess_temp');
    colormap(tile,PIPColorMap);    
    xticks(linspace(0,1,5));
    xticklabels(linspace(0,1,5));
    yticks(linspace(0,1,5));
    yticklabels(linspace(0,1,5));
    set(gca,'PlotBoxAspectRatio',[1 1 1]);
    set(gca,'YDir','normal','FontSize',12,'FontWeight','bold');
    xlabel('p (resident)','FontSize',14,'FontWeight','bold'); 
    ylabel('p (mutant)','FontSize',14,'FontWeight','bold');

end


%% Add horizontal lines on the steady state plots
nexttile(1)
yline(MaximaPoint{1}(1),'LineWidth',1.5,'Color','r','LineStyle','-');
nexttile(3);
yline(MaximaPoint{2}(1),'LineWidth',1.5,'Color','r','LineStyle','-');
nexttile(5);
yline(MaximaPoint{3}(1),'LineWidth',1.5,'Color','r','LineStyle','-');

%% Add Annotations
nexttile(1);
%text(.07235,.91,'$(\gamma_{L_{max}}^{(12hr)},p_{L_{max}}^{(12hr)})$','interpreter','latex','FontSize',12,'FontWeight','bold','HorizontalAlignment','right','Color','k','BackgroundColor','white');
text(.95,.17235,'$(p^*,\gamma^*)$','interpreter','latex','FontSize',12,'FontWeight','bold','HorizontalAlignment','right','Color','k','BackgroundColor','white');
text(MaximaPoint{1}(2),MaximaPoint{1}(1),'$\longrightarrow$','interpreter','latex','FontSize',12,'FontWeight','bold','HorizontalAlignment','right','Color','k','Rotation',-40);
nexttile(2);
text(.25,.65,'+','FontSize',30,'Color',[1 1 1],'FontWeight','bold','HorizontalAlignment','center');
text(.85,.65,'+','FontSize',30,'Color',[1 1 1],'FontWeight','bold','HorizontalAlignment','center');
text(.6,.85,'-','FontSize',30,'Color',0*[1 1 1],'FontWeight','bold','HorizontalAlignment','center');
text(.6,.25,'-','FontSize',30,'Color',0*[1 1 1],'FontWeight','bold','HorizontalAlignment','center');
nexttile(3);
%text(.0505,.82,'$(\gamma_{L_{max}}^{(16hr)},p_{L_{max}}^{(16hr)})$','interpreter','latex','FontSize',12,'FontWeight','bold','HorizontalAlignment','right','Color','k','BackgroundColor','white');
text(.89,.1905,'$(p^*,\gamma^*)$','interpreter','latex','FontSize',12,'FontWeight','bold','HorizontalAlignment','right','Color','k','BackgroundColor','white');
text(MaximaPoint{2}(2),MaximaPoint{2}(1),'$\longrightarrow$','interpreter','latex','FontSize',12,'FontWeight','bold','HorizontalAlignment','right','Color','k','Rotation',-90);
nexttile(4);
text(.15,.5,'+','FontSize',30,'Color',[1 1 1],'FontWeight','bold','HorizontalAlignment','center');
text(.75,.5,'+','FontSize',30,'Color',[1 1 1],'FontWeight','bold','HorizontalAlignment','center');
text(.5,.75,'-','FontSize',30,'Color',0*[1 1 1],'FontWeight','bold','HorizontalAlignment','center');
text(.5,.15,'-','FontSize',30,'Color',0*[1 1 1],'FontWeight','bold','HorizontalAlignment','center');
nexttile(5);
%text(.0505,.52,'$(\gamma_{L_{max}}^{(24hr)},p_{L_{max}}^{(24hr)})$','interpreter','latex','FontSize',12,'FontWeight','bold','HorizontalAlignment','right','Color','k','BackgroundColor','white');
%text(.52,.0505,'$(p^*,\gamma^*)$','interpreter','latex','FontSize',12,'FontWeight','bold','HorizontalAlignment','right','Color','k','BackgroundColor','white');
text(.585,.1905,'$(p^*,\gamma^*)$','interpreter','latex','FontSize',12,'FontWeight','bold','HorizontalAlignment','right','Color','k','BackgroundColor','white');
text(MaximaPoint{3}(2),MaximaPoint{3}(1),'$\longrightarrow$','interpreter','latex','FontSize',12,'FontWeight','bold','HorizontalAlignment','right','Color','k','Rotation',-90);
nexttile(6);
text(.15,.3,'+','FontSize',30,'Color',[1 1 1],'FontWeight','bold','HorizontalAlignment','center');
text(.6,.3,'+','FontSize',30,'Color',[1 1 1],'FontWeight','bold','HorizontalAlignment','center');
text(.3,.75,'-','FontSize',30,'Color',0*[1 1 1],'FontWeight','bold','HorizontalAlignment','center');
text(.3,.15,'-','FontSize',30,'Color',0*[1 1 1],'FontWeight','bold','HorizontalAlignment','center');

%% Add plot labels
nexttile(2);
text(-2.45,1.02,'(A)','FontSize',14,'FontWeight','bold');
nexttile(2);
text(-.25,1.02,'(B)','FontSize',14,'FontWeight','bold');
nexttile(4);
text(-2.45,1.02,'(C)','FontSize',14,'FontWeight','bold');
nexttile(4);
text(-.25,1.02,'(D)','FontSize',14,'FontWeight','bold');
nexttile(6);
text(-2.45,1.02,'(E)','FontSize',14,'FontWeight','bold');
nexttile(6);
text(-.25,1.02,'(F)','FontSize',14,'FontWeight','bold');
   
%% Add plot titles
nexttile(1);
title('Cycle period = 12hr','FontSize',18,'FontWeight','bold','Position',[4.5 1.05]);
nexttile(3);
title('Cycle period = 16hr','FontSize',18,'FontWeight','bold','Position',[4.5 1.05]);
nexttile(5);
title('Cycle period = 24hr','FontSize',18,'FontWeight','bold','Position',[4.5 1.05]);

%% Add labels for the invasion plots

annotation('rectangle',[.855,.9037,.02,.02],'FaceColor',[1 1 1]);
annotation('textbox',[.855,.9027,.02,.02],'String','$\mid$','FitBoxToText','on','EdgeColor','none','FontSize',10,'FontWeight','bold','Color',0*[1 1 1],'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','middle');
annotation('rectangle',[.8796,.9037,.02,.02],'FaceColor',0*[1 1 1]);
annotation('textbox',[.8796,.9027,.02,.02],'String','+','FitBoxToText','on','EdgeColor','none','FontSize',18,'FontWeight','bold','Color',[1 1 1],'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','middle');

annotation('rectangle',[.855,.6086,.02,.02],'FaceColor',[1 1 1]);
annotation('textbox',[.855,.6076,.02,.02],'String','$\mid$','FitBoxToText','on','EdgeColor','none','FontSize',10,'FontWeight','bold','Color',0*[1 1 1],'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','middle');
annotation('rectangle',[.8796,.6086,.02,.02],'FaceColor',0*[1 1 1]);
annotation('textbox',[.8796,.6076,.02,.02],'String','+','FitBoxToText','on','EdgeColor','none','FontSize',18,'FontWeight','bold','Color',[1 1 1],'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','middle');

annotation('rectangle',[.855,.3112,.02,.02],'FaceColor',[1 1 1]);
annotation('textbox',[.855,.3102,.02,.02],'String','$\mid$','FitBoxToText','on','EdgeColor','none','FontSize',10,'FontWeight','bold','Color',0*[1 1 1],'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','middle');
annotation('rectangle',[.8796,.3112,.02,.02],'FaceColor',0*[1 1 1]);
annotation('textbox',[.8796,.3102,.02,.02],'String','+','FitBoxToText','on','EdgeColor','none','FontSize',18,'FontWeight','bold','Color',[1 1 1],'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','middle');

nexttile(2);
text(1.16,.9,'Mutant invades','FontSize',14,'FontWeight','bold','Rotation',-90);
text(1.06,.9,'Mutant invasion fails','FontSize',14,'FontWeight','bold','Rotation',-90);


nexttile(4);
text(1.16,.9,'Mutant invades','FontSize',14,'FontWeight','bold','Rotation',-90);
text(1.06,.9,'Mutant invasion fails','FontSize',14,'FontWeight','bold','Rotation',-90);

nexttile(6);
text(1.16,.9,'Mutant invades','FontSize',14,'FontWeight','bold','Rotation',-90);
text(1.06,.9,'Mutant invasion fails','FontSize',14,'FontWeight','bold','Rotation',-90);



%% Save Figure
filename = dir('..\\Figures\\ChangingCyclePeriod*');

if isempty(filename)
    filename = '..\\Figures\\ChangingCyclePeriod_v1.eps';
else
    filename = [filename(end).folder '\' filename(end).name];
    version = extractBetween(filename,"_v",".");
    version = version{1};
    version = str2num(version);
    filename = [extractBefore(filename,['v' num2str(version)]) 'v' num2str(version+1) '.eps'];
end
%exportgraphics(h,filename,"BackgroundColor",'none','ContentType','vector');

