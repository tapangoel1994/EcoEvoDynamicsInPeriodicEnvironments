%% Code to generate figures S1, S2 and S3

%% Date created: 6/13/2024
%% Author: Tapan Goel

close all; 
clear all;

addpath('utils/');
addpath('lib/');
    colorpalette;
    fixedparameters;

%% simulation parameters:
params.T = 24; % hours
params.t_vals = transpose(0:params.dt:params.T); % time
NCycles = 20;
%% filter parameters

q_L = 0;
q_V = 0;

%% initial conditions
R0 = 1e2; %initial resource amount in ug/mL ( 500 mL flask)
S0 = 1e7; %Initial concentration of susceptibles in flask (per mL)
V01= 1e4; %initial concentration of virus in flask (per mL)
V02 = 0;
x0 = [R0 S0 zeros(1,6) V01 V02];


q_LV = [0 .2;.1 .1; .2 0]; %% each row represents a different filtration condition

figuretitles = {"only virions", "virions and lysogens","only lysogens"};

for conditions = 1:3 %% iterate over different filtration conditions
    
    TransferMatrix = diag([q_R q_S q_E q_E q_I q_I q_LV(conditions,1) q_LV(conditions,1) q_LV(conditions,2) q_LV(conditions,2)]);
    
    
    %% Obligately lytic
    % lysogen probability and induction rate
    params.p = [0 0];
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
    params.p = [1 0];
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
    params.p = [.5 0];
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
    
    
    
    %% Make figure
    
    supfig3_1 = figure('Position',[100 0 1260 800]);
    t_main = tiledlayout(3,5,'TileSpacing','loose','Padding','tight');
    
        
    %%%% Plot lytic strategy

    nexttile(t_main,1,[1 2]);
    semilogy(T_lytic(1:3*length(t_vals)),Y_lytic(1:3*length(t_vals),1),'LineWidth',2,'Color',linecolors.R,'LineStyle',linestyle.lytic); % plot lytic R
    hold on;
    semilogy(T_lytic(1:3*length(t_vals)),Y_lytic(1:3*length(t_vals),2),'LineWidth',2,'Color',linecolors.S,'LineStyle',linestyle.lytic); % plot lytic S
    semilogy(T_lytic(1:3*length(t_vals)),Y_lytic(1:3*length(t_vals),3),'LineWidth',2,'Color',linecolors.E,'LineStyle',linestyle.lytic); % plot lytic E
    semilogy(T_lytic(1:3*length(t_vals)),Y_lytic(1:3*length(t_vals),5),'LineWidth',2,'Color',linecolors.I,'LineStyle',linestyle.lytic); % plot lytic I
    semilogy(T_lytic(1:3*length(t_vals)),Y_lytic(1:3*length(t_vals),7),'LineWidth',2,'Color',linecolors.L,'LineStyle',linestyle.lytic); % plot lytic L
    semilogy(T_lytic(1:3*length(t_vals)),Y_lytic(1:3*length(t_vals),9),'LineWidth',2,'Color',linecolors.V,'LineStyle',linestyle.lytic); % plot lytic V
    
    semilogy(T_lytic(1:length(t_vals):3*length(t_vals)),Y_lytic(1:length(t_vals):3*length(t_vals),1),'LineWidth',1.5,'Marker',marker.lytic,'Color',linecolors.R,'LineStyle','none');
    semilogy(T_lytic(1:length(t_vals):3*length(t_vals)),Y_lytic(1:length(t_vals):3*length(t_vals),2),'LineWidth',1.5,'Marker',marker.lytic,'Color',linecolors.S,'LineStyle','none');
    semilogy(T_lytic(1:length(t_vals):3*length(t_vals)),Y_lytic(1:length(t_vals):3*length(t_vals),3),'LineWidth',1.5,'Marker',marker.lytic,'Color',linecolors.E,'LineStyle','none');
    semilogy(T_lytic(1:length(t_vals):3*length(t_vals)),Y_lytic(1:length(t_vals):3*length(t_vals),5),'LineWidth',1.5,'Marker',marker.lytic,'Color',linecolors.I,'LineStyle','none');
    semilogy(T_lytic(1:length(t_vals):3*length(t_vals)),Y_lytic(1:length(t_vals):3*length(t_vals),7),'LineWidth',1.5,'Marker',marker.lytic,'Color',linecolors.L,'LineStyle','none');
    semilogy(T_lytic(1:length(t_vals):3*length(t_vals)),Y_lytic(1:length(t_vals):3*length(t_vals),9),'LineWidth',1.5,'Marker',marker.lytic,'Color',linecolors.V,'LineStyle','none');

    xline(.99*params.T,'--','LineWidth',1.5,'Color',.8*[1 1 1]);
    xline(1.99*params.T,'--','LineWidth',1.5,'Color',.5*[1 1 1]);
    xline(2.99*params.T,'--','LineWidth',1.5,'Color',.2*[1 1 1]);
    ylim([1e-2 1e10]);
    xlim([0 3.1*params.T]);
    set(gca,'YMinorTick','off','Box','off','XTick',0:8:3*params.T,'XTickLabel',[],'YTick',logspace(-2,10,4));
    set(gca,'FontSize',16,'FontWeight','bold','FontName','Times');
    text(-13,5e10,'(A)','FontSize',16,'FontWeight','bold');
    hold off;
    
    
    nexttile(t_main,3,[1 3]);
    semilogy(0:NCycles-1,Y_lytic(1:length(t_vals):end,1),'LineWidth',1,'Color',linecolors.R,'LineStyle',linestyle.lytic,'Marker',marker.lytic);
    hold on;
    semilogy(0:NCycles-1,Y_lytic(1:length(t_vals):end,2),'LineWidth',1,'Color',linecolors.S,'LineStyle',linestyle.lytic,'Marker',marker.lytic);
    semilogy(0:NCycles-1,Y_lytic(1:length(t_vals):end,3),'LineWidth',1,'Color',linecolors.E,'LineStyle',linestyle.lytic,'Marker',marker.lytic);
    semilogy(0:NCycles-1,Y_lytic(1:length(t_vals):end,5),'LineWidth',1,'Color',linecolors.I,'LineStyle',linestyle.lytic,'Marker',marker.lytic);
    semilogy(0:NCycles-1,Y_lytic(1:length(t_vals):end,7),'LineWidth',1,'Color',linecolors.L,'LineStyle',linestyle.lytic,'Marker',marker.lytic);
    semilogy(0:NCycles-1,Y_lytic(1:length(t_vals):end,9),'LineWidth',1,'Color',linecolors.V,'LineStyle',linestyle.lytic,'Marker',marker.lytic);
    xline(.95,'--','LineWidth',1.5,'Color',.8*[1 1 1]);
    xline(1.95,'--','LineWidth',1.5,'Color',.5*[1 1 1]);
    xline(2.95,'--','LineWidth',1.5,'Color',.2*[1 1 1]);
    ylim([1e-2 1e10]);
    set(gca,'YMinorTick','off','Box','off','XTick',1:2:NCycles,'XTickLabel',[],'YTick',logspace(-2,10,4),'YTickLabel',[]);
    set(gca,'FontSize',16,'FontWeight','bold','FontName','Times');
    
    text(-1,5e10,'(B)','FontSize',16,'FontWeight','bold');
    hold off;
    


    %%% Plot temperate strategy

    nexttile(t_main,6,[1 2]);
    semilogy(T_temperate(1:3*length(t_vals)),Y_temperate(1:3*length(t_vals),1),'LineWidth',2,'Color',linecolors.R); % plot temperate R
    hold on;
    semilogy(T_temperate(1:3*length(t_vals)),Y_temperate(1:3*length(t_vals),2),'LineWidth',2,'Color',linecolors.S); % plot temperate S
    semilogy(T_temperate(1:3*length(t_vals)),Y_temperate(1:3*length(t_vals),3),'LineWidth',2,'Color',linecolors.E); % plot temperate E
    semilogy(T_temperate(1:3*length(t_vals)),Y_temperate(1:3*length(t_vals),5),'LineWidth',2,'Color',linecolors.I); % plot temperate I
    semilogy(T_temperate(1:3*length(t_vals)),Y_temperate(1:3*length(t_vals),7),'LineWidth',2,'Color',linecolors.L); % plot temperate L
    semilogy(T_temperate(1:3*length(t_vals)),Y_temperate(1:3*length(t_vals),9),'LineWidth',2,'Color',linecolors.V); % plot temperate V

    semilogy(T_temperate(1:length(t_vals):3*length(t_vals)),Y_temperate(1:length(t_vals):3*length(t_vals),1),'LineWidth',1.5,'Marker',marker.temperate,'Color',linecolors.R,'LineStyle','none');
    semilogy(T_temperate(1:length(t_vals):3*length(t_vals)),Y_temperate(1:length(t_vals):3*length(t_vals),2),'LineWidth',1.5,'Marker',marker.temperate,'Color',linecolors.S,'LineStyle','none');
    semilogy(T_temperate(1:length(t_vals):3*length(t_vals)),Y_temperate(1:length(t_vals):3*length(t_vals),3),'LineWidth',1.5,'Marker',marker.temperate,'Color',linecolors.E,'LineStyle','none');
    semilogy(T_temperate(1:length(t_vals):3*length(t_vals)),Y_temperate(1:length(t_vals):3*length(t_vals),5),'LineWidth',1.5,'Marker',marker.temperate,'Color',linecolors.I,'LineStyle','none');
    semilogy(T_temperate(1:length(t_vals):3*length(t_vals)),Y_temperate(1:length(t_vals):3*length(t_vals),7),'LineWidth',1.5,'Marker',marker.temperate,'Color',linecolors.L,'LineStyle','none');
    semilogy(T_temperate(1:length(t_vals):3*length(t_vals)),Y_temperate(1:length(t_vals):3*length(t_vals),9),'LineWidth',1.5,'Marker',marker.temperate,'Color',linecolors.V,'LineStyle','none');

    xline(params.T,'--','LineWidth',1.5,'Color',.8*[1 1 1]);
    xline(2*params.T,'--','LineWidth',1.5,'Color',.5*[1 1 1]);
    xline(3*params.T,'--','LineWidth',1.5,'Color',.2*[1 1 1]);
    ylim([1e-2 1e10]);
    xlim([0 3.1*params.T]);
    set(gca,'YMinorTick','off','Box','off','XTick',0:8:3*params.T,'XTickLabel',[],'YTick',logspace(-2,10,4));
    set(gca,'FontSize',16,'FontWeight','bold','FontName','Times');
    text(-13,5e10,'(C)','FontSize',16,'FontWeight','bold');
    hold off;
    
    
    nexttile(t_main,8,[1 3]);
    semilogy(0:NCycles-1,Y_temperate(1:length(t_vals):end,1),'LineWidth',1,'Color',linecolors.R,'Marker',marker.temperate);
    hold on;
    semilogy(0:NCycles-1,Y_temperate(1:length(t_vals):end,2),'LineWidth',1,'Color',linecolors.S,'Marker',marker.temperate);
    semilogy(0:NCycles-1,Y_temperate(1:length(t_vals):end,3),'LineWidth',1,'Color',linecolors.E,'Marker',marker.temperate);
    semilogy(0:NCycles-1,Y_temperate(1:length(t_vals):end,5),'LineWidth',1,'Color',linecolors.I,'Marker',marker.temperate);
    semilogy(0:NCycles-1,Y_temperate(1:length(t_vals):end,7),'LineWidth',1,'Color',linecolors.L,'Marker',marker.temperate);
    semilogy(0:NCycles-1,Y_temperate(1:length(t_vals):end,9),'LineWidth',1,'Color',linecolors.V,'Marker',marker.temperate);
    xline(.95,'--','LineWidth',1.5,'Color',.8*[1 1 1]);
    xline(1.95,'--','LineWidth',1.5,'Color',.5*[1 1 1]);
    xline(2.95,'--','LineWidth',1.5,'Color',.2*[1 1 1]);
    ylim([1e-2 1e10]);
    set(gca,'YMinorTick','off','Box','off','XTick',1:2:NCycles,'XTickLabel',[],'YTick',logspace(-2,10,4),'YTickLabel',[]);
    set(gca,'FontSize',16,'FontWeight','bold','FontName','Times');
    text(-1,5e10,'(D)','FontSize',16,'FontWeight','bold');
    hold off;
    
    
    %%% Plot lysogenic strategy
   
    nexttile(t_main,11,[1 2]);
    semilogy(T_lysogenic(1:3*length(t_vals)),Y_lysogenic(1:3*length(t_vals),1),'LineWidth',2,'Color',linecolors.R); % plot lysogenic R
    hold on;
    semilogy(T_lysogenic(1:3*length(t_vals)),Y_lysogenic(1:3*length(t_vals),2),'LineWidth',2,'Color',linecolors.S); % plot lysogenic S
    semilogy(T_lysogenic(1:3*length(t_vals)),Y_lysogenic(1:3*length(t_vals),3),'LineWidth',2,'Color',linecolors.E); % plot lysogenic E
    semilogy(T_lysogenic(1:3*length(t_vals)),Y_lysogenic(1:3*length(t_vals),5),'LineWidth',2,'Color',linecolors.I); % plot lysogenic I
    semilogy(T_lysogenic(1:3*length(t_vals)),Y_lysogenic(1:3*length(t_vals),7),'LineWidth',2,'Color',linecolors.L); % plot lysogenic L
    semilogy(T_lysogenic(1:3*length(t_vals)),Y_lysogenic(1:3*length(t_vals),9),'LineWidth',2,'Color',linecolors.V); % plot lysogenic V

    semilogy(T_lysogenic(1:length(t_vals):3*length(t_vals)),Y_lysogenic(1:length(t_vals):3*length(t_vals),1),'LineWidth',1.5,'Marker',marker.lysogenic,'Color',linecolors.R,'LineStyle','none');
    semilogy(T_lysogenic(1:length(t_vals):3*length(t_vals)),Y_lysogenic(1:length(t_vals):3*length(t_vals),2),'LineWidth',1.5,'Marker',marker.lysogenic,'Color',linecolors.S,'LineStyle','none');
    semilogy(T_lysogenic(1:length(t_vals):3*length(t_vals)),Y_lysogenic(1:length(t_vals):3*length(t_vals),3),'LineWidth',1.5,'Marker',marker.lysogenic,'Color',linecolors.E,'LineStyle','none');
    semilogy(T_lysogenic(1:length(t_vals):3*length(t_vals)),Y_lysogenic(1:length(t_vals):3*length(t_vals),5),'LineWidth',1.5,'Marker',marker.lysogenic,'Color',linecolors.I,'LineStyle','none');
    semilogy(T_lysogenic(1:length(t_vals):3*length(t_vals)),Y_lysogenic(1:length(t_vals):3*length(t_vals),7),'LineWidth',1.5,'Marker',marker.lysogenic,'Color',linecolors.L,'LineStyle','none');
    semilogy(T_lysogenic(1:length(t_vals):3*length(t_vals)),Y_lysogenic(1:length(t_vals):3*length(t_vals),9),'LineWidth',1.5,'Marker',marker.lysogenic,'Color',linecolors.V,'LineStyle','none');

    xline(params.T,'--','LineWidth',1.5,'Color',.8*[1 1 1]);
    xline(2*params.T,'--','LineWidth',1.5,'Color',.5*[1 1 1]);
    xline(3*params.T,'--','LineWidth',1.5,'Color',.2*[1 1 1]);
    ylim([1e-2 1e10]);
    xlim([0 3.1*params.T]);
    set(gca,'YMinorTick','off','Box','off','XTick',0:8:3*params.T,'XTickLabel',0:8:3*params.T,'YTick',logspace(-2,10,4));
    set(gca,'FontSize',16,'FontWeight','bold','FontName','Times');
    xlabel('Time (hr)','FontSize',18,'FontWeight','bold','FontName','Times');
    text(-13,5e10,'(E)','FontSize',16,'FontWeight','bold');
    hold off;
    
    
    nexttile(t_main,13,[1 3]);
    semilogy(0:NCycles-1,Y_lysogenic(1:length(t_vals):end,1),'LineWidth',1,'Color',linecolors.R,'Marker',marker.lysogenic);
    hold on;
    semilogy(0:NCycles-1,Y_lysogenic(1:length(t_vals):end,2),'LineWidth',1,'Color',linecolors.S,'Marker',marker.lysogenic);
    semilogy(0:NCycles-1,Y_lysogenic(1:length(t_vals):end,3),'LineWidth',1,'Color',linecolors.E,'Marker',marker.lysogenic);
    semilogy(0:NCycles-1,Y_lysogenic(1:length(t_vals):end,5),'LineWidth',1,'Color',linecolors.I,'Marker',marker.lysogenic);
    semilogy(0:NCycles-1,Y_lysogenic(1:length(t_vals):end,7),'LineWidth',1,'Color',linecolors.L,'Marker',marker.lysogenic);
    semilogy(0:NCycles-1,Y_lysogenic(1:length(t_vals):end,9),'LineWidth',1,'Color',linecolors.V,'Marker',marker.lysogenic);
    xline(.95,'--','LineWidth',1.5,'Color',.8*[1 1 1]);
    xline(1.95,'--','LineWidth',1.5,'Color',.5*[1 1 1]);
    xline(2.95,'--','LineWidth',1.5,'Color',.2*[1 1 1]);
    ylim([1e-2 1e10]);
    set(gca,'YMinorTick','off','Box','off','XTick',1:2:NCycles,'XTickLabel',2:2:NCycles,'YTick',logspace(-2,10,4),'YTickLabel',[]);
    set(gca,'FontSize',16,'FontWeight','bold','FontName','Times');
    xlabel('Number of cycles','FontSize',18,'FontWeight','bold','FontName','Times');
    text(-1,5e10,'(F)','FontSize',16,'FontWeight','bold');
    hold off;
    
    
    nexttile(6)
    ylabel('Density (mL$^{-1}$)','FontSize',18,'FontWeight','bold','FontName','Times');
    
     %% Figure title
    titlehandle = title(t_main,{sprintf("Filtration condition: %s pass through",figuretitles{conditions});""},'FontSize',20,'FontWeight','bold','FontName','Times','interpreter','latex');
    
    %% Add legends
    nexttile(t_main,1,[1 2]);
    legend('R','S','E','I','L','V','FontName','Times','FontWeight','bold','FontSize',14,'Position',[0.0897,0.9154,0.2977,0.03299],'orientation','horizontal');
    legend('Box','off');
    nexttile(t_main,3,[1 3]);
    hold on;
    lh(1) = plot(nan,nan,'dk');
    lh(2) = plot(nan,nan,'sk');
    lh(3) = plot(nan,nan,'ok');
    L{1} = 'Obligately lytic';
    L{2} = 'Temperate';
    L{3} = 'Obligately lysogenic';
    legend(lh,L,'Position',[0.5114,0.9154,0.4297,0.03379],'Orientation','horizontal','Box','off');
    legend('Box','off');

   
    %% Save Figure
    filename = dir(['../Figures/FigureS' num2str(conditions) '*']);
    
    if isempty(filename)
        filename = ['../Figures/FigureS' num2str(conditions) '_v1.eps'];
    else
        filename = ['../Figures/' filename(end).name];
        version = extractBetween(filename,"_v",".");
        version = version{1};
        version = str2num(version);
        filename = [extractBefore(filename,['v' num2str(version)]) 'v' num2str(version+1) '.eps'];
    end
    filename1 = [filename(1:end-4) '.png'];

    exportgraphics(supfig3_1,filename,"BackgroundColor",'none','ContentType','vector');
    exportgraphics(supfig3_1,filename1,"BackgroundColor",'white');
end