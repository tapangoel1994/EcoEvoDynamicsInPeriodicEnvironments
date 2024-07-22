%% Code to simulate and plot single cycle for 3 different viral strategies over 24
%% hours. This script generates Figure 2 in the manuscript.

%% Date created: 2/23/2024
%% Modified: 6/25/2024
%% Modified: 6/26/2024
%% Author: Tapan Goel

close all; 
clear all;

%% load fixed parameters and color scheme
addpath('utils\');
addpath('lib\');
    colorpalette;
    fixedparameters;


%% simulation parameters:
params.T = 24; % hours
params.t_vals = transpose(0:params.dt:params.T); % time


%% initial conditions
R0 = 1e2; %initial resource amount in ug/mL
S0 = 1e7; %Initial concentration of susceptibles in flask (per mL)
V01= 1e4; %initial concentration of virus in flask (per mL)
V02 = 0;
x0 = [R0 S0 zeros(1,6) V01 V02];

PGamma = [0 0;.5 .083;1 0]; %% each row is the (p,gamma) for a particular strategy
StrategyLabels = {'Obligately lytic','Temperate','Obligately lysogenic'};

h = figure('Position',[100 0 400 1200]);
t = tiledlayout(4,1,"TileSpacing",'compact','Padding','tight');

for i = 1:length(PGamma)

    params.p = [PGamma(i,1) 0];
    params.gamma = [PGamma(i,2) 0];

    [t_vals,y] = ode113(@ODE_RSEILV_2Species,params.t_vals,x0,options,params);

    nexttile(i);
    semilogy(t_vals,y(:,1),'LineWidth',3,'Color',linecolors.R);  %Plot R
    hold on;
    semilogy(t_vals,y(:,2),'LineWidth',3,'Color',linecolors.S); %Plot S
    semilogy(t_vals,y(:,3),'LineWidth',3,'Color',linecolors.E); %Plot E1
    semilogy(t_vals,y(:,5),'LineWidth',3,'Color',linecolors.I); %Plot I1
    semilogy(t_vals,y(:,7),'LineWidth',3,'Color',linecolors.L); %Plot L1
    semilogy(t_vals,y(:,9),'LineWidth',3,'Color',linecolors.V); %Plot V1
    ylim([1e-2 1e10]);
    xlim([0 max(t_vals)]);
    set(gca,'YMinorTick','off','Box','off','XTick',0:4:24,'XTickLabel',[],'YTick',logspace(-2,10,5));
    text(11,1e9,StrategyLabels{i},'FontSize',14,'FontWeight','bold');
    set(gca,'FontSize',14,'FontWeight','bold');
   
    if i == 1
        lgd1 = legend('R','S','E','I','L','V');
        lgd1.Position = [0.5527    0.7710    0.1590    0.1062];
        legend('boxoff');
    end

    %if i == 2
        ylabel('Density $($mL$^{-1})$','FontSize',16,'FontWeight','bold');
    %end
    nexttile(4);
    semilogy(t_vals,sum(y(:,3:10),2),'LineWidth',3,'Color','k','LineStyle',linestyle.(Strategies{i}));
    hold on;

end

nexttile(4);
ylim([1e2 1e8]);
xlim([0 max(t_vals)]);
set(gca,'YMinorTick','off','Box','off','XTick',0:4:24,'XTickLabel',0:4:24,'YTick',logspace(2,8,4));
lgd = legend("Obligately lytic","Temperate","Obligately lysogenic",'Location','best');
legend('boxoff');
lgd.Position = [0.1927    0.2120    0.4964    0.0546];
set(gca,'FontSize',14,'FontWeight','bold');
ylabel('Total viral genomes (mL$^{-1}$)','FontSize',16,'FontWeight','bold','Position',[-3.5169 3.87e4]);
xlabel('Time (hr)','FontSize',16,'FontWeight','bold');

%% Add panel labels
nexttile(1);
text(-5.5,10^10,'(A)','FontSize',16,'FontWeight','bold');
nexttile(2);
text(-5.5,10^10,'(B)','FontSize',16,'FontWeight','bold');
nexttile(3);
text(-5.5,10^10,'(C)','FontSize',16,'FontWeight','bold');
nexttile(4);
text(-5.5,10^8,'(D)','FontSize',16,'FontWeight','bold');
%% Save Figure
filename = dir('..\Figure\Figure2*');

if isempty(filename)
    filename = '..\Figures\Figure2_v1.eps';
else
    filename = ['..\Figures\' filename(end).name];
    version = extractBetween(filename,"_v",".");
    version = version{1};
    version = str2num(version);
    filename = [extractBefore(filename,['v' num2str(version)]) 'v' num2str(version+1) '.eps'];    
end
filename1 = [filename(1:end-4) '.png'];

exportgraphics(h,filename,"BackgroundColor",'none','ContentType','vector');
exportgraphics(h,filename1,"BackgroundColor",'white');
