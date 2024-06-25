%%Code to simulate and plot single cycle for 3 different viral strategies over 24
%%hours

%% Date created: 2/23/2024
%% Modified: 6/25/2024
%% Author: Tapan Goel

%close all; 
%clear all;

addpath('utils\');
load('lib\colorpalette.mat');
load('lib\fixedparameters.mat');


%% simulation parameters:
params.T = 24; % hours
params.t_vals = transpose(0:params.dt:params.T); % time

%% filter parameters
p_R = 1;
p_S = 0;
p_E = 0;
p_I = 0;
p_L = 0.2;
p_V = 0.0;
TransferMatrix = diag([p_R p_S p_E p_E p_I p_I p_L p_L p_V p_V]);


%% initial conditions
R0 = 1e2; %initial resource amount in ug/mL ( 500 mL flask)
S0 = 1e7; %Initial concentration of susceptibles in flask (per mL)
V01= 1e4; %initial concentration of virus in flask (per mL)
V02 = 0;
x0 = [R0 S0 zeros(1,6) V01 V02];


%% Obligate lytic
% lysogen probability and induction rate
params.q = [0 0];
params.gamma = [0 0];
[t_vals, y_lytic] = ode113(@ODE_RSEILV_2Species, params.t_vals, x0, options, params);


%% Purely lysogenic
params.q = [1 0];
params.gamma = [0 0];
[t_vals, y_lysogenic] = ode113(@ODE_RSEILV_2Species, params.t_vals, x0, options, params);


%% Temperate
params.q = [.5 0];
params.gamma = [.083 0];
[t_vals, y_temperate] = ode113(@ODE_RSEILV_2Species, params.t_vals, x0, options, params);


%% Plot all three time series

h = figure('Position',[100 0 3*600 1200]);
t = tiledlayout(2,6,"TileSpacing",'tight');

nexttile(1,[1 2]);
    semilogy(t_vals,y_lytic(:,1),'LineWidth',3,'Color',linecolors.R);  %Plot R
    hold on;
    semilogy(t_vals,y_lytic(:,2),'LineWidth',3,'Color',linecolors.S); %Plot S
    semilogy(t_vals,y_lytic(:,3),'LineWidth',3,'Color',linecolors.E); %Plot E1
    semilogy(t_vals,y_lytic(:,5),'LineWidth',3,'Color',linecolors.I); %Plot I1
    semilogy(t_vals,y_lytic(:,7),'LineWidth',3,'Color',linecolors.L); %Plot L1
    semilogy(t_vals,y_lytic(:,9),'LineWidth',3,'Color',linecolors.V); %Plot V1
    ylim([1e-2 1e10]);
    xlim([0 max(t_vals)]);
    set(gca,'YMinorTick','off','Box','off','XTick',0:4:24,'XTickLabel',0:4:24,'YTick',logspace(-2,10,7));
    title('Obligately lytic');
    set(gca,'FontSize',22,'FontWeight','bold','FontName','Times');
    ylabel('Density (mL^{-1})','FontSize',26,'FontWeight','bold','FontName','Times');
    legend('R','S','E','I','L','V','location','best');

nexttile(3,[1 2]);
    semilogy(t_vals,y_temperate(:,1),'LineWidth',3,'Color',linecolors.R);  %Plot R
    hold on;
    semilogy(t_vals,y_temperate(:,2),'LineWidth',3,'Color',linecolors.S); %Plot S         
    semilogy(t_vals,y_temperate(:,3),'LineWidth',3,'Color',linecolors.E); %Plot E1
    semilogy(t_vals,y_temperate(:,5),'LineWidth',3,'Color',linecolors.I); %Plot I1
    semilogy(t_vals,y_temperate(:,7),'LineWidth',3,'Color',linecolors.L); %Plot L1
    semilogy(t_vals,y_temperate(:,9),'LineWidth',3,'Color',linecolors.V); %Plot V1
    ylim([1e-2 1e10]);
    xlim([0 max(t_vals)]);
    set(gca,'YMinorTick','off','Box','off','XTick',0:4:24,'XTickLabel',0:4:24,'YTick',logspace(-2,10,7),'YTickLabel',[]);
    title('Temperate');
    set(gca,'FontSize',22,'FontWeight','bold','FontName','Times');
    xlabel(t,'Time (hr)','FontSize',26,'FontWeight','bold','FontName','Times');

nexttile(5,[1 2]);
    semilogy(t_vals,y_lysogenic(:,1),'LineWidth',3,'Color',linecolors.R);  %Plot R
    hold on;
    semilogy(t_vals,y_lysogenic(:,2),'LineWidth',3,'Color',linecolors.S); %Plot S
    semilogy(t_vals,y_lysogenic(:,3),'LineWidth',3,'Color',linecolors.E); %Plot E1
    semilogy(t_vals,y_lysogenic(:,5),'LineWidth',3,'Color',linecolors.I); %Plot I1
    semilogy(t_vals,y_lysogenic(:,7),'LineWidth',3,'Color',linecolors.L); %Plot L1
    semilogy(t_vals,y_lysogenic(:,9),'LineWidth',3,'Color',linecolors.V); %Plot V1
    ylim([1e-2 1e10]);
    xlim([0 max(t_vals)]);
    set(gca,'YMinorTick','off','Box','off','XTick',0:4:24,'XTickLabel',0:4:24,'YTick',logspace(-2,10,7),'YTickLabel',[]);
    title('Obligately lysogenic');
    set(gca,'FontSize',22,'FontWeight','bold','FontName','Times');

nexttile(8,[1 4])
    semilogy(t_vals,sum(y_lytic(:,3:10),2),'LineWidth',3,'Color','k','LineStyle',linestyle.lytic);
    hold on;
    semilogy(t_vals,sum(y_temperate(:,3:10),2),'LineWidth',3,'Color','k','LineStyle',linestyle.temperate);
    semilogy(t_vals,sum(y_lysogenic(:,3:10),2),'LineWidth',3,'Color','k','LineStyle',linestyle.lysogenic);
    ylim([1e-2 1e10]);
    xlim([0 max(t_vals)]);
    set(gca,'YMinorTick','off','Box','off','XTick',0:4:24,'XTickLabel',0:4:24,'YTick',logspace(-2,10,7));
    set(gca,'FontSize',22,'FontWeight','bold','FontName','Times');
    ylabel('Total viral genomes (mL^{-1})','FontSize',26,'FontWeight','bold','FontName','Times');
    legend("Obligately lytic","Temperate","Obligately lysogenic")

 

 
