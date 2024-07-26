%% Script generates steady state viral genome density plots as function of Q
%% for the value of gamma over which the invasion dynamics have been calculated.
%% this is to check that the excluded regions in the invasion analysis make sense. 

%% This script generates figure S5.


%% Author: Tapan Goel
%% Date: 7/18/2024

clear all;
close all;

addpath("lib/");
addpath("utils/");
addpath("../Data/");
colorpalette;
fixedparameters;


%% Generate genome density plots corresponding to figure 4

CyclePeriod = 24;
Gamma = logspace(-3,0,51);
P = linspace(0,1,51);
NumNodes = 12;

%%Add missing simulation parameters
params.T = CyclePeriod; % hours
params.t_vals = transpose(0:params.dt:params.T); % time

q_LV = [0 .2;.2 0];
S0 = 1e7;
Va_0 = 1e4;
Vb_0 = 10*criticaldensitythreshold;

for index = 1:2
    if isfile(sprintf("../Data/SteadyState_CyclePeriod=%.1f,S0=%1.e,V0=%1.e,q_L=%.1f,q_V=%.1f.mat",CyclePeriod,S0,Va_0,q_LV(index,1),q_LV(index,2)))
        load(sprintf("SteadyState_CyclePeriod=%.1f,S0=%1.e,V0=%1.e,q_L=%.1f,q_V=%.1f.mat",CyclePeriod,S0,Va_0,q_LV(index,1),q_LV(index,2)));
        SteadyState{index} = SteadyStateDensity;
        CyclesSteadyState{index} = SSCycles;
    else
        [SteadyState{index}, CyclesSteadyState{index}] = PopulationSteadyStateFunction(CyclePeriod,q_LV(index,1),q_LV(index,2),Gamma,P,NumNodes,1,params);
    end
end

% Plot the steady state viral genome densities of resident along the fixed
% gamma line

h = figure('Renderer','painters','Position',[760 200 1225 951/2]);
t = tiledlayout(1,2,'Padding','loose','TileSpacing','loose');

inset_position = [0.1600    0.3000    0.1994    0.4000; .62 .33 .1994 .4];
colorlabel_position = [2.9109 4.4107e+06 0; 2.9109 3e2 0];
for index = 1:2
    
    SteadyState_temp = SteadyState{index};
    [~,i] = min(abs(Gamma - 0.0832));
    
    nexttile(t,index);
    semilogy(P,sum(SteadyState_temp(:,i,3:2:10),3),'LineWidth',2,'Color','k','Marker','o','MarkerFaceColor','k','MarkerSize',4);
    hold on;
    yline(100*criticaldensitythreshold,'Color',[.5 .5 .5],'LineWidth',1.5,'LineStyle','--');
    ylim([1e-2 10*max(sum(SteadyState_temp(:,i,3:2:10),3))]);

    set(gca,'YMinorTick','off','Box','off','XTick',0:.25:1,'YTick',logspace(-2,10,7),'XColor','r','LineWidth',2);
    set(gca,'FontSize',16,'FontWeight','bold','FontName','Times');
    hold off;

    %% Add inset figure
    ax_inset{index} = axes('Position',inset_position(index,:));
    imagesc(P,Gamma,sum(SteadyState_temp(:,:,3:2:10),3)');    
    hold on;
    %plot(Gamma(i),Q(j),'*k','MarkerSize',5,'LineWidth',2);
    contour(P,Gamma,sum(SteadyState_temp(:,:,3:2:10),3)','k');
    yline(.0832,'LineWidth',1.5,'Color','r');
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
    xlabel( '$p$: Integration probability');
    pbaspect(ax_inset{index},[1 1 1]);
    hold off;

end

xlabel(t,'$p$: Integration probability','FontSize',18,'Interpreter','latex');
ylabel(t,{'Viral genome density (mL$^{-1}$)'},'FontSize',18,'Interpreter','latex');
nexttile(t,1);
title('Filtrate: virions only','FontSize',20);
nexttile(t,2);
title('Filtrate: lysogens only','FontSize',20);

%%Add panel labels
%% Add plot labels
nexttile(t,1);
text(-.21,3e10,'(A)','FontSize',20,'FontWeight','bold');
text(1.15,3e10,'(B)','FontSize',20,'FontWeight','bold');


%% Save Figure
filename = dir('../Figures/FigureS5*');

if isempty(filename)
    filename = '../Figures/FigureS5_v1.eps';
else
    filename = ['../Figures/' filename(end).name];
    version = extractBetween(filename,"_v",".");
    version = version{1};
    version = str2num(version);
    filename = [extractBefore(filename,['v' num2str(version)]) 'v' num2str(version+1) '.eps'];
end
filename1 = [filename(1:end-4) '.png'];

exportgraphics(h,filename,"BackgroundColor",'none','ContentType','vector');
exportgraphics(h,filename1,"BackgroundColor",'white');