%%Code to generate figure 4 of the manuscript: Steady State Density and PIP
%%plots for passaging only virions and only lysogens.

%% Date created: 05/29/2024
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

p_LV = [.2 0;.1 .1;0 .2];
S0 = 1e7;
V01 = 1e4;


for index = 1:3
    if isfile(sprintf("..\\Data\\CyclePeriod=%.1f,S0=%1.e,V0=%1.e,p_L=%.1f,p_V=%.1f.mat",CyclePeriod,S0,V01,p_LV(index,1),p_LV(index,2)))
        load(sprintf("..\\Data\\CyclePeriod=%.1f,S0=%1.e,V0=%1.e,p_L=%.1f,p_V=%.1f.mat",CyclePeriod,S0,V01,p_LV(index,1),p_LV(index,2)));
        SteadyState{index} = SteadyStateDensity;
        CyclesSteadyState{index} = SSCycles;
    else
        [SteadyState{index}, CyclesSteadyState{index}] = PopulationSteadyStateFunction(CyclePeriod,p_LV(index,1),p_LV(index,2),Gamma,Q,NumNodes,1,params);
    end
end

for index = 1:3
    
    SteadyState_temp = SteadyState{index};
    [M,I] = max(SteadyState_temp(:,:,7),[],"all","linear");
    [j,i]=ind2sub([length(Q),length(Gamma)],I);
    
    
    InvasionVariable = [Q' .0832*ones(size(Q'))];
    %InvasionVariable = [Q' Gamma(i)*ones(size(Q'))];
    if isfile(sprintf("..\\Data\\InvasionCyclePeriod=%.1f,S0=%1.e,V0=%1.e,p_L=%.1f,p_V=%.1f.mat",CyclePeriod,S0,V01,p_LV(index,1),p_LV(index,2)))
        load(sprintf("..\\Data\\InvasionCyclePeriod=%.1f,S0=%1.e,V0=%1.e,p_L=%.1f,p_V=%.1f.mat",CyclePeriod,S0,V01,p_LV(index,1),p_LV(index,2)));
        Invasion{index} = InvasionSteadyStateDensity;
        CyclesInvasion{index} = InvasionSSCycles;
    else
        [Invasion{index}, CyclesInvasion{index}] = InvasionSteadyStateFunction(CyclePeriod,p_LV(index,1),p_LV(index,2),InvasionVariable,NumNodes,1,params);
    end


end

%% Plot 
h = figure('Renderer','painters','Position',[760 200 1225 951]);
t = tiledlayout(2,2,'Padding','loose','TileSpacing','loose');

for index = 3:-2:1
    
    SteadyState_temp = SteadyState{index};
    SteadyState_temp = squeeze(sum(SteadyState_temp(:,:,3:10),3));
    M = max(SteadyState_temp,[],"all","linear");
    [r,c] = find(SteadyState_temp == M);
    j = max(r);
    i = max(c);

    Invasion_temp = Invasion{index};
    resident = squeeze(sum(Invasion_temp(:,:,3:2:end),3));
    mutant = squeeze(sum(Invasion_temp(:,:,4:2:end),3));

    reswin = (resident > 1e-1/params.flask_volume) & (mutant < 1e-1/params.flask_volume);
    mutwin = (resident < 1e-1/params.flask_volume) & (mutant > 1e-1/params.flask_volume);
    coexist = (resident > 1e-1/params.flask_volume) & (mutant > 1e-1/params.flask_volume);
    alldie = (resident < 1e-1/params.flask_volume) & (mutant < 1e-1/params.flask_volume);

    InvasionSuccess = zeros(length(Q),length(Q));
    %InvasionSuccess = mutwin*2 + coexist*1 - reswin*1 - alldie*2;
    InvasionSuccess = double(mutant > resident) - double(mutant < resident);
    InvasionSuccess(1:length(Q)+1:end) = 0;
    
    %% Plot steady state densities
    tile = nexttile;
    imagesc(Gamma,Q,SteadyState_temp);
    hold on;
    plot(Gamma(i),Q(j),'*r','MarkerSize',5,'LineWidth',2);
    
    
    
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
    
    %yticklabels({'0','.5','1'});
    set(gca,'XScale','log','YScale','linear','ColorScale','log','YDir','normal','FontWeight','bold','FontSize',14); 
    xlabel('$\gamma$: Inducation rate (hr$^{-1}$)','FontSize',16,"FontWeight",'bold');
    ylabel( '$p$: Integration Probability','FontSize',16,"FontWeight",'bold');
    hold off;

    %% Plot PIP
    tile = nexttile;
    imagesc(Q,Q,InvasionSuccess');
    colormap(tile,flipud(gray));    
    xticks(linspace(0,1,5));
    xticklabels(linspace(0,1,5));
    yticks(linspace(0,1,5));
    yticklabels(linspace(0,1,5));
    set(gca,'YDir','normal','FontSize',14,'FontWeight','bold');
    
    xlabel('p (resident)','FontSize',16,'FontWeight','bold'); 
    ylabel('p (mutant)','FontSize',16,'FontWeight','bold');

end


%% Add vertical lines on the steady state plots
nexttile(1)
xline(.0832,'LineWidth',1.5,'Color','k','LineStyle','--');
nexttile(3)
xline(.0832,'LineWidth',1.5,'Color','k','LineStyle','--');
%% Add Annotations
nexttile(1);
text(0.6210,0.0508,'$(\gamma_m,p_m)$','interpreter','latex','FontSize',14,'FontWeight','bold','HorizontalAlignment','right','Color','red','BackgroundColor','white');
text(1,0,'$\longrightarrow$','interpreter','latex','FontSize',14,'FontWeight','bold','HorizontalAlignment','right','Color','red','Rotation',-15);
nexttile(2);
text(.75,.25,'+','FontSize',30,'Color',[1 1 1],'FontWeight','bold','HorizontalAlignment','center');
text(.25,.75,'-','FontSize',30,'Color',0*[1 1 1],'FontWeight','bold','HorizontalAlignment','center');
nexttile(3);
text(.0495,.5175,'$(\gamma_m,p_m)$','interpreter','latex','FontSize',14,'FontWeight','bold','HorizontalAlignment','right','Color','red','BackgroundColor','white');
text(.0832,.52,'$\longrightarrow$','interpreter','latex','FontSize',14,'FontWeight','bold','HorizontalAlignment','right','Color','red');
nexttile(4);
text(.15,.3,'+','FontSize',30,'Color',[1 1 1],'FontWeight','bold','HorizontalAlignment','center');
text(.6,.3,'+','FontSize',30,'Color',[1 1 1],'FontWeight','bold','HorizontalAlignment','center');
text(.3,.75,'-','FontSize',30,'Color',0*[1 1 1],'FontWeight','bold','HorizontalAlignment','center');
text(.3,.15,'-','FontSize',30,'Color',0*[1 1 1],'FontWeight','bold','HorizontalAlignment','center');

%% Add plot labels
nexttile(1);
text(10^-3.6,1,'(A)','FontSize',16,'FontWeight','bold');
nexttile(2);
text(-.18,1,'(B)','FontSize',16,'FontWeight','bold');
nexttile(3);
text(10^-3.6,1,'(C)','FontSize',16,'FontWeight','bold');
nexttile(4);
text(-.18,1,'(D)','FontSize',16,'FontWeight','bold');
   
%% Add plot titles
nexttile(1);
title('Filtrate: virions only','FontSize',20,'FontWeight','bold','Position',[4.5 1.05]);
nexttile(3);
title('Filtrate: lysogens only','FontSize',20,'FontWeight','bold','Position',[4.5 1.05]);

%% Add labels for the invasion plots

annotation('rectangle',[.9124,.9037,.02,.02],'FaceColor',[1 1 1]);
annotation('textbox',[.9124,.9027,.02,.02],'String','$\mid$','FitBoxToText','on','EdgeColor','none','FontSize',12,'FontWeight','bold','Color',0*[1 1 1],'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','middle');
annotation('rectangle',[.942,.9037,.02,.02],'FaceColor',0*[1 1 1]);
annotation('textbox',[.942,.9027,.02,.02],'String','+','FitBoxToText','on','EdgeColor','none','FontSize',20,'FontWeight','bold','Color',[1 1 1],'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','middle');

annotation('rectangle',[.9124,.428,.02,.02],'FaceColor',[1 1 1]);
annotation('textbox',[.9124,.427,.02,.02],'String','$\mid$','FitBoxToText','on','EdgeColor','none','FontSize',12,'FontWeight','bold','Color',0*[1 1 1],'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','middle');
annotation('rectangle',[.942,.428,.02,.02],'FaceColor',0*[1 1 1]);
annotation('textbox',[.942,.427,.02,.02],'String','+','FitBoxToText','on','EdgeColor','none','FontSize',20,'FontWeight','bold','Color',[1 1 1],'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','middle');

nexttile(2);
text(1.16,.92,'Mutant invades','FontSize',16,'FontWeight','bold','Rotation',-90);
text(1.06,.92,'Mutant invasion fails','FontSize',16,'FontWeight','bold','Rotation',-90);


nexttile(4);
text(1.16,.92,'Mutant invades','FontSize',16,'FontWeight','bold','Rotation',-90);
text(1.06,.92,'Mutant invasion fails','FontSize',16,'FontWeight','bold','Rotation',-90);


%% Save Figure
filename = dir('..\\Figures\\SteadyStatePlusInvasion*');

if isempty(filename)
    filename = '..\Figures\SteadyStatePlusInvasion_v1.eps';
else
    filename = filename(end).name;
    version = extractBetween(filename,"_v",".");
    version = version{1};
    version = str2num(version);
    filename = [extractBefore(filename,num2str(version)) num2str(version+1) '.eps'];
end
%saveas(h,filename,'epsc');

