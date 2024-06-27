%%Code to simulate and plot multiple cycle for 3 different viral strategies over 24
%%hours

%% Date created: 6/13/2024
%% Author: Tapan Goel

close all; 
clear all;

addpath('Utils\');
load('lib\colorpalette.mat');
load('lib\fixedparameters.mat');


%% simulation parameters:
params.T = 24; % hours
params.t_vals = transpose(0:params.dt:params.T); % time
NCycles = 20;
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


%% Initiate figure
h_main = figure('Position',[100 0 1260 800]);
t_main = tiledlayout(3,5,'TileSpacing','compact','Padding','compact');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Only virions passage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_L = 0.0;
p_V = 0.2;
TransferMatrix = diag([p_R p_S p_E p_E p_I p_I p_L p_L p_V p_V]);


%% Obligately lytic
% lysogen probability and induction rate
params.q = [0 0];
params.gamma = [0 0];
T_lytic = [];
Y_lytic = [];
x0 = [R0 S0 zeros(1,6) V01 V02];
for iter = 1:NCycles
    [t_vals, y_lytic] = ode113(@ODE_RSEILV_2Species, params.t_vals, x0, options, params);
    T_lytic = [T_lytic; params.T*(iter-1) + t_vals];
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
    T_lysogenic = [T_lysogenic; params.T*(iter-1) + t_vals];
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
    T_temperate = [T_temperate; params.T*(iter-1) + t_vals];
    Y_temperate = [Y_temperate;y_temperate];
    x0 = [R0 S0 zeros(1,8)] + y_temperate(end,:)*TransferMatrix;
end



figure(h_main);
nexttile(t_main,1,[1 2]);
semilogy(T_lytic(1:3*length(t_vals)),sum(Y_lytic(1:3*length(t_vals),3:end),2),'LineWidth',1.5,'Color','k','LineStyle',linestyle.lytic);
hold on;
semilogy(T_temperate(1:3*length(t_vals)),sum(Y_temperate(1:3*length(t_vals),3:end),2),'LineWidth',1.5,'Color','k','LineStyle',linestyle.temperate);
semilogy(T_lysogenic(1:3*length(t_vals)),sum(Y_lysogenic(1:3*length(t_vals),3:end),2),'LineWidth',1.5,'Color','k','LineStyle',linestyle.lysogenic);

semilogy(T_lytic(1:length(t_vals):3*length(t_vals)),sum(Y_lytic(1:length(t_vals):3*length(t_vals),3:end),2),'Marker',marker.lytic,'Color','k','LineStyle','none');
semilogy(T_temperate(1:length(t_vals):3*length(t_vals)),sum(Y_temperate(1:length(t_vals):3*length(t_vals),3:end),2),'Marker',marker.temperate,'Color','k','LineStyle','none');
semilogy(T_lysogenic(1:length(t_vals):3*length(t_vals)),sum(Y_lysogenic(1:length(t_vals):3*length(t_vals),3:end),2),'Marker',marker.lysogenic,'Color','k','LineStyle','none');

xline(0.99*params.T,'--','LineWidth',1.5,'Color',.8*[1 1 1]);
xline(1.99*params.T,'--','LineWidth',1.5,'Color',.5*[1 1 1]);
xline(2.99*params.T,'--','LineWidth',1.5,'Color',.2*[1 1 1]);
ylim([1e-2 1e10]);
xlim([0 3.1*params.T]);
set(gca,'YMinorTick','off','Box','off','XTick',0:8:3*params.T,'XTickLabel',[],'YTick',logspace(-2,10,4));
set(gca,'FontSize',16,'FontWeight','bold','FontName','Times');
text(-13,1e10,'(A)','FontSize',16,'FontWeight','bold');
text(.3*params.T,1e0,'Cycle 1','FontSize',16,'FontWeight','bold');
text(1.3*params.T,1e0,'Cycle 2','FontSize',16,'FontWeight','bold');
text(2.3*params.T,1e0,'Cycle 3','FontSize',16,'FontWeight','bold');
hold off;


nexttile(t_main,3,[1 3]);
semilogy(0:NCycles-1,sum(Y_lytic(1:length(t_vals):end,3:end),2),'LineWidth',1,'Color','k','LineStyle',linestyle.lytic,'Marker',marker.lytic);
hold on;
semilogy(0:NCycles-1,sum(Y_temperate(1:length(t_vals):end,3:end),2),'LineWidth',1,'Color','k','LineStyle',linestyle.temperate,'Marker',marker.temperate);
semilogy(0:NCycles-1,sum(Y_lysogenic(1:length(t_vals):end,3:end),2),'LineWidth',1,'Color','k','LineStyle',linestyle.lysogenic,'Marker',marker.lysogenic);
xline(.95,'--','LineWidth',1.5,'Color',.8*[1 1 1]);
xline(1.95,'--','LineWidth',1.5,'Color',.5*[1 1 1]);
xline(2.95,'--','LineWidth',1.5,'Color',.2*[1 1 1]);

ylim([1e-2 1e10]);
set(gca,'YMinorTick','off','Box','off','XTick',1:2:NCycles,'XTickLabel',[],'YTick',logspace(-2,10,4),'YTickLabel',[]);
set(gca,'FontSize',16,'FontWeight','bold','FontName','Times');
legend('Obligately lytic','Temperate','Obligately lysogenic','FontName','Times','FontWeight','bold','FontSize',14,'location','best');
legend('boxoff');
text(-1,1e10,'(B)','FontSize',16,'FontWeight','bold');
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Lysogens and virions passage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_L = 0.1;
p_V = 0.1;
TransferMatrix = diag([p_R p_S p_E p_E p_I p_I p_L p_L p_V p_V]);


%% Obligately lytic
% lysogen probability and induction rate
params.q = [0 0];
params.gamma = [0 0];
T_lytic = [];
Y_lytic = [];
x0 = [R0 S0 zeros(1,6) V01 V02];
for iter = 1:NCycles
    [t_vals, y_lytic] = ode113(@ODE_RSEILV_2Species, params.t_vals, x0, options, params);
    T_lytic = [T_lytic; params.T*(iter-1) + t_vals];
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
    T_lysogenic = [T_lysogenic; params.T*(iter-1) + t_vals];
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
    T_temperate = [T_temperate; params.T*(iter-1) + t_vals];
    Y_temperate = [Y_temperate;y_temperate];
    x0 = [R0 S0 zeros(1,8)] + y_temperate(end,:)*TransferMatrix;
end


figure(h_main);
nexttile(t_main,6,[1 2]);
semilogy(T_lytic(1:3*length(t_vals)),sum(Y_lytic(1:3*length(t_vals),3:end),2),'LineWidth',1.5,'Color','k','LineStyle',linestyle.lytic);
hold on;
semilogy(T_temperate(1:3*length(t_vals)),sum(Y_temperate(1:3*length(t_vals),3:end),2),'LineWidth',1.5,'Color','k','LineStyle',linestyle.temperate);
semilogy(T_lysogenic(1:3*length(t_vals)),sum(Y_lysogenic(1:3*length(t_vals),3:end),2),'LineWidth',1.5,'Color','k','LineStyle',linestyle.lysogenic);

semilogy(T_lytic(1:length(t_vals):3*length(t_vals)),sum(Y_lytic(1:length(t_vals):3*length(t_vals),3:end),2),'Marker',marker.lytic,'Color','k','LineStyle','none');
semilogy(T_temperate(1:length(t_vals):3*length(t_vals)),sum(Y_temperate(1:length(t_vals):3*length(t_vals),3:end),2),'Marker',marker.temperate,'Color','k','LineStyle','none');
semilogy(T_lysogenic(1:length(t_vals):3*length(t_vals)),sum(Y_lysogenic(1:length(t_vals):3*length(t_vals),3:end),2),'Marker',marker.lysogenic,'Color','k','LineStyle','none');

xline(0.99*params.T,'--','LineWidth',1.5,'Color',.8*[1 1 1]);
xline(1.99*params.T,'--','LineWidth',1.5,'Color',.5*[1 1 1]);
xline(2.99*params.T,'--','LineWidth',1.5,'Color',.2*[1 1 1]);
ylim([1e-2 1e10]);
xlim([0 3.1*params.T]);
set(gca,'YMinorTick','off','Box','off','XTick',0:8:3*params.T,'XTickLabel',[],'YTick',logspace(-2,10,4));
set(gca,'FontSize',16,'FontWeight','bold','FontName','Times');
text(-13,1e10,'(C)','FontSize',16,'FontWeight','bold');
hold off;


nexttile(t_main,8,[1 3]);
semilogy(0:NCycles-1,sum(Y_lytic(1:length(t_vals):end,3:end),2),'LineWidth',1,'Color','k','LineStyle',linestyle.lytic,'Marker',marker.lytic);
hold on;
semilogy(0:NCycles-1,sum(Y_temperate(1:length(t_vals):end,3:end),2),'LineWidth',1,'Color','k','LineStyle',linestyle.temperate,'Marker',marker.temperate);
semilogy(0:NCycles-1,sum(Y_lysogenic(1:length(t_vals):end,3:end),2),'LineWidth',1,'Color','k','LineStyle',linestyle.lysogenic,'Marker',marker.lysogenic);
xline(.95,'--','LineWidth',1.5,'Color',.8*[1 1 1]);
xline(1.95,'--','LineWidth',1.5,'Color',.5*[1 1 1]);
xline(2.95,'--','LineWidth',1.5,'Color',.2*[1 1 1]);
ylim([1e-2 1e10]);
set(gca,'YMinorTick','off','Box','off','XTick',1:2:NCycles,'XTickLabel',[],'YTick',logspace(-2,10,4),'YTickLabel',[]);
set(gca,'FontSize',16,'FontWeight','bold','FontName','Times');
text(-1,1e10,'(D)','FontSize',16,'FontWeight','bold');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Only lysogens passage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_L = 0.2;
p_V = 0;
TransferMatrix = diag([p_R p_S p_E p_E p_I p_I p_L p_L p_V p_V]);


%% Obligately lytic
% lysogen probability and induction rate
params.q = [0 0];
params.gamma = [0 0];
T_lytic = [];
Y_lytic = [];
x0 = [R0 S0 zeros(1,6) V01 V02];
for iter = 1:NCycles
    [t_vals, y_lytic] = ode113(@ODE_RSEILV_2Species, params.t_vals, x0, options, params);
    T_lytic = [T_lytic; params.T*(iter-1) + t_vals];
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
    T_lysogenic = [T_lysogenic; params.T*(iter-1) + t_vals];
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
    T_temperate = [T_temperate; params.T*(iter-1) + t_vals];
    Y_temperate = [Y_temperate;y_temperate];
    x0 = [R0 S0 zeros(1,8)] + y_temperate(end,:)*TransferMatrix;
end

figure(h_main);
nexttile(t_main,11,[1 2]);
semilogy(T_lytic(1:3*length(t_vals)),sum(Y_lytic(1:3*length(t_vals),3:end),2),'LineWidth',1.5,'Color','k','LineStyle',linestyle.lytic);
hold on;
semilogy(T_temperate(1:3*length(t_vals)),sum(Y_temperate(1:3*length(t_vals),3:end),2),'LineWidth',1.5,'Color','k','LineStyle',linestyle.temperate);
semilogy(T_lysogenic(1:3*length(t_vals)),sum(Y_lysogenic(1:3*length(t_vals),3:end),2),'LineWidth',1.5,'Color','k','LineStyle',linestyle.lysogenic);

semilogy(T_lytic(1:length(t_vals):3*length(t_vals)),sum(Y_lytic(1:length(t_vals):3*length(t_vals),3:end),2),'Marker',marker.lytic,'Color','k','LineStyle','none');
semilogy(T_temperate(1:length(t_vals):3*length(t_vals)),sum(Y_temperate(1:length(t_vals):3*length(t_vals),3:end),2),'Marker',marker.temperate,'Color','k','LineStyle','none');
semilogy(T_lysogenic(1:length(t_vals):3*length(t_vals)),sum(Y_lysogenic(1:length(t_vals):3*length(t_vals),3:end),2),'Marker',marker.lysogenic,'Color','k','LineStyle','none');

xline(0.99*params.T,'--','LineWidth',1.5,'Color',.8*[1 1 1]);
xline(1.99*params.T,'--','LineWidth',1.5,'Color',.5*[1 1 1]);
xline(2.99*params.T,'--','LineWidth',1.5,'Color',.2*[1 1 1]);
ylim([1e-2 1e10]);
xlim([0 3.1*params.T]);
set(gca,'YMinorTick','off','Box','off','XTick',0:8:3*params.T,'XTickLabel',0:8:3*params.T,'YTick',logspace(-2,10,4));
set(gca,'FontSize',16,'FontWeight','bold','FontName','Times');
xlabel('Times (hr)','FontSize',18,'FontWeight','bold','FontName','Times');
text(-13,1e10,'(E)','FontSize',16,'FontWeight','bold');
hold off;


nexttile(t_main,13,[1 3]);
semilogy(0:NCycles-1,sum(Y_lytic(1:length(t_vals):end,3:end),2),'LineWidth',1,'Color','k','LineStyle',linestyle.lytic,'Marker',marker.lytic);
hold on;
semilogy(0:NCycles-1,sum(Y_temperate(1:length(t_vals):end,3:end),2),'LineWidth',1,'Color','k','LineStyle',linestyle.temperate,'Marker',marker.temperate);
semilogy(0:NCycles-1,sum(Y_lysogenic(1:length(t_vals):end,3:end),2),'LineWidth',1,'Color','k','LineStyle',linestyle.lysogenic,'Marker',marker.lysogenic);
xline(.95,'--','LineWidth',1.5,'Color',.8*[1 1 1]);
xline(1.95,'--','LineWidth',1.5,'Color',.5*[1 1 1]);
xline(2.95,'--','LineWidth',1.5,'Color',.2*[1 1 1]);
ylim([1e-2 1e10]);
set(gca,'YMinorTick','off','Box','off','XTick',1:2:NCycles,'XTickLabel',2:2:NCycles,'YTick',logspace(-2,10,4),'YTickLabel',[]);
set(gca,'FontSize',16,'FontWeight','bold','FontName','Times');
xlabel('Number of cycles','FontSize',18,'FontWeight','bold','FontName','Times');
text(-1,1e10,'(F)','FontSize',16,'FontWeight','bold');
hold off;

nexttile(6);
ylabel('Total viral genome density (mL$^{-1})$','FontSize',18,'FontWeight','bold','FontName','Times');