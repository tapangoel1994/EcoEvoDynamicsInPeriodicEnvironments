%%Code to simulate and plot multiple cycle for 3 different viral strategies over 24
%%hours

%% Date created: 3/14/2024
%% Author: Tapan Goel

close all; 
clear all;

addpath('Utils');

%% Life history parameters (units of hours, micrograms and mL). 
params.J = 0; %ug/mL-h
params.conversion_efficiency = 5e-7; %ug/cell
params.d_R = 0; % per hour
params.mu_max = 1.2; % per hour
params.R_in = 4; %ug/mL
params.alpha_l = 0;
params.alpha_e = 0;
params.alpha_i = 0;

params.d_S = .2; %per hour
params.d_E = .2; %per hour
params.d_L = .2; %per hour
params.d_I = .2; %per hour
params.m = 1/24; %per hour

params.phi = 3.4e-10; %mL/hr
params.lambda = 2; %per hour
params.eta = 1; %per hour
params.bet = 50;




%% simulation parameters:
params.flask_volume = 500; %volume in mL
params.dt = 1/30; % hours
params.T = 24; % hours
params.t_vals = transpose(0:params.dt:params.T); % time
NCycles = 5;

%% filter parameters
p_S = 0;
p_E = 0;
p_I = 0;
p_L = 0.2;
p_V = 0;
TransferMatrix = diag([0 p_S p_E p_E p_I p_I p_L p_L p_V p_V]);

%% Numerical method related parameters
options = odeset('AbsTol',1e-8,'RelTol',1e-8,'NonNegative',1:10); %Options for the ODE function call
steadystatethresh = 1e-1/params.flask_volume; % concentration difference below which two concentrations are treated as identical


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
T_lytic = [];
Y_lytic = [];
x0 = [R0 S0 zeros(1,6) V01 V02];
for iter = 1:NCycles
    [t_vals, y_lytic] = ode113(@ODE_RSEILV_2Species, params.t_vals, x0, options, params);
    T_lytic = [T_lytic; 24*(iter-1) + t_vals];
    Y_lytic = [Y_lytic;y_lytic];
    x0 = [R0 S0 zeros(1,8)] + y_lytic(end,:)*TransferMatrix;
end

%% Purely lysogenic
params.q = [1 0];
params.gamma = [0 0];
T_lysogenic = [];
Y_lysogenic = [];
x0 = [R0 S0 zeros(1,6) V01 V02];
for iter = 1:NCycles
    [t_vals, y_lysogenic] = ode113(@ODE_RSEILV_2Species, params.t_vals, x0, options, params);
    T_lysogenic = [T_lysogenic; 24*(iter-1) + t_vals];
    Y_lysogenic = [Y_lysogenic;y_lysogenic];
    x0 = [R0 S0 zeros(1,8)] + y_lysogenic(end,:)*TransferMatrix;
end


%% Temperate
params.q = [.5 0];
params.gamma = [.083 0];
x0 = [R0 S0 zeros(1,6) V01 V02];
T_temperate = [];
Y_temperate = [];
for iter = 1:NCycles
    [t_vals, y_temperate] = ode113(@ODE_RSEILV_2Species, params.t_vals, x0, options, params);
    T_temperate = [T_temperate; 24*(iter-1) + t_vals];
    Y_temperate = [Y_temperate;y_temperate];
    x0 = [R0 S0 zeros(1,8)] + y_temperate(end,:)*TransferMatrix;
end


%% Plot all three time series
markerstyle = {'+','*','.','x','s'};
    linecolors.S = '#000000';
    linecolors.E = '#ffa600';
    linecolors.I = '#a74ac7';
    linecolors.L = '#ff0000';
    linecolors.V = '#00cbff';

h = figure('Position',[100 0 1200 1200]);
t = tiledlayout(3,1);

nexttile(1);
   
    semilogy(T_lytic,Y_lytic(:,1),'LineWidth',3,'Color',[.5 .5 .5]);  %Plot R
    hold on;
    semilogy(T_lytic,Y_lytic(:,2),'LineWidth',3,'Color',linecolors.S); %Plot S
    
    semilogy(T_lytic,Y_lytic(:,3),'LineWidth',3,'Color',linecolors.E); %Plot E1
    semilogy(T_lytic,Y_lytic(:,5),'LineWidth',3,'Color',linecolors.I); %Plot I1
    semilogy(T_lytic,Y_lytic(:,7),'LineWidth',3,'Color',linecolors.L); %Plot L1
    semilogy(T_lytic,Y_lytic(:,9),'LineWidth',3,'Color',linecolors.V); %Plot V1
     xline(params.T:params.T:params.T*(NCycles-1),'--','Color',[.75 .75 .75],'LineWidth',2);
    ylim([1e-2 1e10]);
    xlim([0 max(T_lytic)]);
    set(gca,'YMinorTick','off','Box','off','XTickLabel',0:params.T:NCycles*params.T,'XTick',0:params.T:NCycles*params.T,'YTick',logspace(-2,10,7));
    title('Obligate lytic');
    legend('R','S','E','I','L','V','location','best');
    set(gca,'FontSize',18);
    set(gca,"FontWeight",'bold');   
    ylabel(t,'Density (mL^{-1})','FontSize',22,'FontWeight','bold');

nexttile(2);
    semilogy(T_temperate,Y_temperate(:,1),'LineWidth',3,'Color',[.5 .5 .5]);  %Plot R
    hold on;
    semilogy(T_temperate,Y_temperate(:,2),'LineWidth',3,'Color',linecolors.S); %Plot S
                        
    semilogy(T_temperate,Y_temperate(:,3),'LineWidth',3,'Color',linecolors.E); %Plot E1
    semilogy(T_temperate,Y_temperate(:,5),'LineWidth',3,'Color',linecolors.I); %Plot I1
    semilogy(T_temperate,Y_temperate(:,7),'LineWidth',3,'Color',linecolors.L); %Plot L1
    semilogy(T_temperate,Y_temperate(:,9),'LineWidth',3,'Color',linecolors.V); %Plot V1
     xline(params.T:params.T:params.T*(NCycles-1),'--','Color',[.75 .75 .75],'LineWidth',2);
    ylim([1e-2 1e10]);
    xlim([0 max(T_lytic)]);
    set(gca,'YMinorTick','off','Box','off','XTickLabel',0:params.T:NCycles*params.T,'XTick',0:params.T:NCycles*params.T,'YTick',logspace(-2,10,7));
    title('Temperate');
    set(gca,'FontSize',18);
    set(gca,"FontWeight",'bold');  
    xlabel(t,'Time (hr)','FontSize',22,'FontWeight','bold');

nexttile(3);
    semilogy(T_lysogenic,Y_lysogenic(:,1),'LineWidth',3,'Color',[.5 .5 .5]);  %Plot R
    hold on;
    semilogy(T_lysogenic,Y_lysogenic(:,2),'LineWidth',3,'Color',linecolors.S); %Plot S
                        
    semilogy(T_lysogenic,Y_lysogenic(:,3),'LineWidth',3,'Color',linecolors.E); %Plot E1
    semilogy(T_lysogenic,Y_lysogenic(:,5),'LineWidth',3,'Color',linecolors.I); %Plot I1
    semilogy(T_lysogenic,Y_lysogenic(:,7),'LineWidth',3,'Color',linecolors.L); %Plot L1
    semilogy(T_lysogenic,Y_lysogenic(:,9),'LineWidth',3,'Color',linecolors.V); %Plot V1
     xline(params.T:params.T:params.T*(NCycles-1),'--','Color',[.75 .75 .75],'LineWidth',2);
    ylim([1e-2 1e10]);
    xlim([0 max(T_lytic)]);
    set(gca,'YMinorTick','off','Box','off','XTickLabel',0:params.T:NCycles*params.T,'XTick',0:params.T:NCycles*params.T,'YTick',logspace(-2,10,7));
    title('Lysogenic');
    set(gca,'FontSize',18);
    set(gca,"FontWeight",'bold');  