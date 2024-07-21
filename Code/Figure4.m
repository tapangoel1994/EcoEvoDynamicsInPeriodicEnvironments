%%Code to generate figure 4 of the manuscript: Steady State Density and PIP
%%plots for passaging only virions and only lysogens.

%% Date created: 07/11/2024
%% Author: Tapan Goel

%%Code pulled from scripts for steady state analysis and PIP analysis using
%%parallelization
clear all;
close all;
addpath('utils\');
addpath('lib\');
colorpalette;
fixedparameters;

%% Define variables
CyclePeriod = 24;
Gamma = logspace(-3,0,51);
Q = linspace(0,1,51);
NumNodes = 12;

%% Add missing simulation parameters
params.T = CyclePeriod; % hours
params.t_vals = transpose(0:params.dt:params.T); % time

p_LV = [.2 0;0 .2];
S0 = 1e7;
Va_0 = 1e4;
Vb_0 = 10*criticaldensitythreshold;

for index = 1:2
    if isfile(sprintf("..\\Data\\SteadyState_CyclePeriod=%.1f,S0=%1.e,V0=%1.e,p_L=%.1f,p_V=%.1f.mat",CyclePeriod,S0,Va_0,p_LV(index,1),p_LV(index,2)))
        load(sprintf("..\\Data\\SteadyState_CyclePeriod=%.1f,S0=%1.e,V0=%1.e,p_L=%.1f,p_V=%.1f.mat",CyclePeriod,S0,Va_0,p_LV(index,1),p_LV(index,2)));
        SteadyState{index} = SteadyStateDensity;
        CyclesSteadyState{index} = SSCycles;
    else
        [SteadyState{index}, CyclesSteadyState{index}] = PopulationSteadyStateFunction(CyclePeriod,p_LV(index,1),p_LV(index,2),Gamma,Q,NumNodes,1,params);
    end
end

for index = 1:2
    
    SteadyState_temp = SteadyState{index};
    [M,I] = max(sum(SteadyState_temp(:,:,3:2:10),3),[],"all","linear");
    [j,i]=ind2sub([length(Q),length(Gamma)],I);
    
    
    InvasionVariable = [Q' .0832*ones(size(Q'))];
    
    if isfile(sprintf("..\\Data\\Invasion_CyclePeriod=%.1f,S0=%1.e,V0=%1.e,p_L=%.1f,p_V=%.1f.mat",CyclePeriod,S0,Va_0,p_LV(index,1),p_LV(index,2)))
        load(sprintf("..\\Data\\Invasion_CyclePeriod=%.1f,S0=%1.e,V0=%1.e,p_L=%.1f,p_V=%.1f.mat",CyclePeriod,S0,Va_0,p_LV(index,1),p_LV(index,2)));
        Invasion{index} = InvasionDensity;
        CyclesInvasion{index} = CyclesToInvasion;
        InvasionSuccessMatrix{index} = InvasionMatrix;

    else
        [Invasion{index}, InvasionSuccessMatrix{index}, CyclesInvasion{index}] = InvasionDynamics(CyclePeriod,p_LV(index,1),p_LV(index,2),InvasionVariable,NumNodes,1,params);
        
    end


end

%% Plot 
h = figure('Renderer','painters','Position',[760 200 1225 951]);
t = tiledlayout(2,2,'Padding','loose','TileSpacing','loose');

for index = 2:-1:1
    
    SteadyState_temp = SteadyState{index};
    SteadyState_temp = squeeze(sum(SteadyState_temp(:,:,3:10),3));
    InvasionSuccess_temp = InvasionSuccessMatrix{index};
    M = max(SteadyState_temp,[],"all","linear");
    [r,c] = find(SteadyState_temp == M);
    j = max(r);
    i = max(c);
    
    %% Plot steady state densities
    tile = nexttile;
    imagesc(Gamma,Q,SteadyState_temp);
    hold on;
    plot(Gamma(i),Q(j),'*k','MarkerSize',5,'LineWidth',2);
    
    
    
    contour(Gamma,Q,SteadyState_temp,'k');
    xticks([1e-3 1e-2 1e-1 1e0]);
    xticklabels({'10$^{-3}$','10$^{-2}$','10$^{-1}$','1'});
    yticks(linspace(0,1,5));
    c = colorbar;
    c.Label.String = "Viral genome density (mL$^{-1}$)";
    c.Label.Interpreter = 'latex';
    c.Label.Rotation = -90;
    c.FontSize = 14;
    c.Label.FontSize = 16;
    

    set(gca,'XScale','log','YScale','linear','ColorScale','log','YDir','normal','FontWeight','bold','FontSize',14); 
    xlabel('$\gamma$: Inducation rate (hr$^{-1}$)','FontSize',16,"FontWeight",'bold');
    ylabel( '$p$: Integration Probability','FontSize',16,"FontWeight",'bold');

    c.Label.Position(1) = 4.7;

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
    set(gca,'YDir','normal','FontSize',14,'FontWeight','bold');
    xlabel('p (resident)','FontSize',16,'FontWeight','bold'); 
    ylabel('p (mutant)','FontSize',16,'FontWeight','bold');

end


%% Add vertical lines on the steady state plots
nexttile(1)
xline(.0832,'LineWidth',1.5,'Color','r','LineStyle','-');
nexttile(3)
xline(.0832,'LineWidth',1.5,'Color','r','LineStyle','-');
%% Add Annotations
nexttile(1);
%text(0.6210,0.0508,'$(\gamma_{V_{max}},p_{V_{max}})$','interpreter','latex','FontSize',14,'FontWeight','bold','HorizontalAlignment','right','Color','k','BackgroundColor','white');
text(0.6210,0.0508,'$(\gamma^*,p^*)$','interpreter','latex','FontSize',14,'FontWeight','bold','HorizontalAlignment','right','Color','k','BackgroundColor','white');
text(1,0,'$\longrightarrow$','interpreter','latex','FontSize',14,'FontWeight','bold','HorizontalAlignment','right','Color','k','Rotation',-15);
nexttile(2);
text(.75,.25,'+','FontSize',30,'Color',[1 1 1],'FontWeight','bold','HorizontalAlignment','center');
text(.25,.75,'-','FontSize',30,'Color',0*[1 1 1],'FontWeight','bold','HorizontalAlignment','center');
nexttile(3);
text(.0495,.5175,'$(\gamma^*,p^*)$','interpreter','latex','FontSize',14,'FontWeight','bold','HorizontalAlignment','right','Color','k','BackgroundColor','white');
text(.0832,.52,'$\longrightarrow$','interpreter','latex','FontSize',14,'FontWeight','bold','HorizontalAlignment','right','Color','k');
nexttile(4);
text(.15,.3,'+','FontSize',30,'Color',[1 1 1],'FontWeight','bold','HorizontalAlignment','center');
text(.6,.3,'+','FontSize',30,'Color',[1 1 1],'FontWeight','bold','HorizontalAlignment','center');
text(.3,.75,'-','FontSize',30,'Color',0*[1 1 1],'FontWeight','bold','HorizontalAlignment','center');
text(.3,.15,'-','FontSize',30,'Color',0*[1 1 1],'FontWeight','bold','HorizontalAlignment','center');

%% Add plot labels
nexttile(1);
text(10^-3.6,1,'(A)','FontSize',16,'FontWeight','bold');
nexttile(2);
text(-.25,1,'(B)','FontSize',16,'FontWeight','bold');
nexttile(3);
text(10^-3.6,1,'(C)','FontSize',16,'FontWeight','bold');
nexttile(4);
text(-.25,1,'(D)','FontSize',16,'FontWeight','bold');
   
%% Add plot titles
nexttile(1);
title('Filtrate: virions only','FontSize',20,'FontWeight','bold','Position',[4.5 1.05]);
nexttile(3);
title('Filtrate: lysogens only','FontSize',20,'FontWeight','bold','Position',[4.5 1.05]);

%% Add labels for the invasion plots

annotation('rectangle',[.873,.9037,.02,.02],'FaceColor',[1 1 1]);
annotation('textbox',[.873,.9027,.02,.02],'String','$\mid$','FitBoxToText','on','EdgeColor','none','FontSize',12,'FontWeight','bold','Color',0*[1 1 1],'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','middle');
annotation('rectangle',[.8967,.9037,.02,.02],'FaceColor',0*[1 1 1]);
annotation('textbox',[.8967,.9027,.02,.02],'String','+','FitBoxToText','on','EdgeColor','none','FontSize',20,'FontWeight','bold','Color',[1 1 1],'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','middle');

annotation('rectangle',[.873,.428,.02,.02],'FaceColor',[1 1 1]);
annotation('textbox',[.873,.427,.02,.02],'String','$\mid$','FitBoxToText','on','EdgeColor','none','FontSize',12,'FontWeight','bold','Color',0*[1 1 1],'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','middle');
annotation('rectangle',[.8967,.428,.02,.02],'FaceColor',0*[1 1 1]);
annotation('textbox',[.8967,.427,.02,.02],'String','+','FitBoxToText','on','EdgeColor','none','FontSize',20,'FontWeight','bold','Color',[1 1 1],'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','middle');

nexttile(2);
text(1.18,.92,'Mutant invades','FontSize',16,'FontWeight','bold','Rotation',-90);
text(1.08,.92,'Mutant invasion fails','FontSize',16,'FontWeight','bold','Rotation',-90);


nexttile(4);
text(1.18,.92,'Mutant invades','FontSize',16,'FontWeight','bold','Rotation',-90);
text(1.08,.92,'Mutant invasion fails','FontSize',16,'FontWeight','bold','Rotation',-90);


%% Save Figure
filename = dir('..\Figures\SteadyStatePlusInvasion*');

if isempty(filename)
    filename = '..\Figures\SteadyStatePlusInvasion_v1.eps';
else
    filename = ['..\Figures\' filename(end).name];
    version = extractBetween(filename,"_v",".");
    version = version{1};
    version = str2num(version);
    filename = [extractBefore(filename,['v' num2str(version)]) 'v' num2str(version+1) '.eps'];
end
exportgraphics(h,filename,"BackgroundColor",'none','ContentType','vector');

