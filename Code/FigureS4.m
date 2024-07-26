%%Code to generate time series till invasion steady state for a two virus-host
%%system for a given (q,gamma) pair.

%% this code generates figure S4

%% Date created: 3/22/22024
%% Author: Tapan Goel

close all; 
clear all;

%addpath('..\..\SerialPassage2\');
addpath('utils\');
addpath('lib\');
colorpalette;
fixedparameters;


% lysogen probability and induction rate for the resident and the mutant
% strains.
params.p = [.5 .3];
params.gamma = [.083 .083];


%% simulation parameters:
params.T = 24; % hours
params.t_vals = transpose(0:params.dt:params.T); % time
MaxCycles = 10000;

%% filter parameters
q_S = 0;
q_E = 0;
q_I = 0;
q_L = 0.2;
q_V = 0.0;
TransferMatrix = diag([0 q_S q_E q_E q_I q_I q_L q_L q_V q_V]);

%% Numerical method related parameters
options = odeset('AbsTol',1e-8,'RelTol',1e-8,'NonNegative',1:10); %Options for the ODE function call
steadystatethresh = 1e-1/params.flask_volume; % concentration difference below which two concentrations are treated as identical


%% initial conditions
R0 = 1e2; %initial resource amount in ug/mL ( 500 mL flask)
S0 = 1e7; %Initial concentration of susceptibles in flask (per mL)
V01= 1e4; %initial concentration of virus in flask (per mL)
V02 = 1e2;


%%%% Let resident reach steady state


%% First cycle
x0 = [R0 S0 zeros(1,6) V01 0];
[t_vals, y] = ode113(@ODE_RSEILV_2Species, params.t_vals, x0, options, params);

TimeSeries = y;
%PlotTimeSeries(t_vals,y)

%% Cycles till steady state or till MaxCycles

iter = 1;
steadyrep = 0;

x0 = [R0 S0 zeros(1,8)] + y(end,:)*TransferMatrix;

%% loop continues till:
%1. The initial population for the new epoch becomes (almost) equal to the
%initial population for the previous epoch for atleast 10 cycles OR,
%2. The total number of cycles hits the max number of cycles.

while ( (sum(abs(x0 - y(1,:)) > steadystatethresh) >= 1) || steadyrep < 10) && iter < MaxCycles+10
    
    if((sum(abs(x0 - y(1,:)) > steadystatethresh) < 1))
                steadyrep = steadyrep+1;
    end

iter = iter+1;

[t_vals, y] = ode113(@ODE_RSEILV_2Species, params.t_vals, x0, options, params);
TimeSeries = [TimeSeries;y];
% PlotTimeSeries_2species(t_vals,y);
% drawnow;
% pause(0.1);
x0 = [R0 S0 zeros(1,8)] + y(end,:)*TransferMatrix;
end
Iter = iter;

%%%% Add mutant and reach steady state

%% First cycle
x0 = [R0 S0 zeros(1,6) 0 V02] + y(end,:)*TransferMatrix;
[t_vals, y] = ode113(@ODE_RSEILV_2Species, params.t_vals, x0, options, params);

TimeSeries = [TimeSeries;y];
%PlotTimeSeries(t_vals,y)

%% Cycles till steady state or till MaxCycles

iter = 1;
steadyrep = 0;

x0 = [R0 S0 zeros(1,8)] + y(end,:)*TransferMatrix;

%% loop continues till:
%1. The initial population for the new epoch becomes (almost) equal to the
%initial population for the previous epoch for atleast 10 cycles OR,
%2. The total number of cycles hits the max number of cycles.

while ( (sum(abs(x0 - y(1,:)) > steadystatethresh) >= 1) || steadyrep < 10) && iter < MaxCycles+10
    
    if((sum(abs(x0 - y(1,:)) > steadystatethresh) < 1))
                steadyrep = steadyrep+1;
    end

iter = iter+1;

[t_vals, y] = ode113(@ODE_RSEILV_2Species, params.t_vals, x0, options, params);
TimeSeries = [TimeSeries;y];
% PlotTimeSeries_2species(t_vals,y);
% drawnow;
% pause(0.1);
x0 = [R0 S0 zeros(1,8)] + y(end,:)*TransferMatrix;
end
Iter = Iter+iter;

%% Plot cycle-to-cycle dynamics



    fig = figure('Renderer','painters','Position',[696 155 793.92 1122.24]);
    %fig = figure('WindowState','maximized');
    t = tiledlayout(3,2,'TileSpacing','loose','Padding','loose');
    
    nexttile(1,[1 2]);
    
    
    semilogy(1:4:Iter,TimeSeries(length(t_vals):4*length(t_vals):end,2),'LineWidth',2,'Color',linecolors.S); %Plot S
    hold on;
    semilogy(1:4:Iter,TimeSeries(length(t_vals):4*length(t_vals):end,3),'LineWidth',2,'Color',linecolors.E); %Plot E1
    semilogy(1:4:Iter,TimeSeries(length(t_vals):4*length(t_vals):end,5),'LineWidth',2,'Color',linecolors.I); %Plot I1
    semilogy(1:4:Iter,TimeSeries(length(t_vals):4*length(t_vals):end,7),'LineWidth',2,'Color',linecolors.L); %Plot L1
    semilogy(1:4:Iter,TimeSeries(length(t_vals):4*length(t_vals):end,9),'LineWidth',2,'Color',linecolors.V); %Plot V1

    semilogy(1:4:Iter,TimeSeries(length(t_vals):4*length(t_vals):end,4),'LineWidth',2,'Color',linecolors.E,'LineStyle','--'); %Plot E2  
    semilogy(1:4:Iter,TimeSeries(length(t_vals):4*length(t_vals):end,6),'LineWidth',2,'Color',linecolors.I,'LineStyle','--'); %Plot I2    
    semilogy(1:4:Iter,TimeSeries(length(t_vals):4*length(t_vals):end,8),'LineWidth',2,'Color',linecolors.L,'LineStyle','--'); %Plot L2    
    semilogy(1:4:Iter,TimeSeries(length(t_vals):4*length(t_vals):end,10),'LineWidth',2,'Color',linecolors.V,'LineStyle','--'); %Plot V2


    legend('S','E$_r$','I$_r$','L$_r$','V$_r$','E$_m$','I$_m$','L$_m$','V$_m$','location','best','NumColumns',2);
    ylim([1e-2 1e10]);
   
    
    set(gca,'YMinorTick','off','Box','off');
    set(gca,'YTick',[1e0 1e2 1e4 1e6 1e8 1e10]);
    %title('Total population');
    set(gca,'FontSize',14,'FontWeight','bold','FontName','Times New Roman');
    xlabel('Number of cycles','FontSize',14,"FontWeight",'bold','FontName','Times New Roman');
    ylabel('Density (mL$^{-1}$)','FontSize',14,"FontWeight",'bold','FontName','Times New Roman');

   %%Pick cycles of interest and plot;
   CycleOfInterest = [1,49,101,201];
   
   for i = 1:length(CycleOfInterest)
       nexttile(2+i);

       semilogy(t_vals,TimeSeries(CycleOfInterest(i)*length(t_vals)+1:(CycleOfInterest(i)+1)*length(t_vals),1),'LineWidth',3,'Color',[.5 .5 .5]);  %Plot R
        hold on;
        semilogy(t_vals,TimeSeries(CycleOfInterest(i)*length(t_vals)+1:(CycleOfInterest(i)+1)*length(t_vals),2),'LineWidth',3,'Color',linecolors.S); %Plot S
        
        semilogy(t_vals,TimeSeries(CycleOfInterest(i)*length(t_vals)+1:(CycleOfInterest(i)+1)*length(t_vals),3),'LineWidth',3,'Color',linecolors.E); %Plot E1
        semilogy(t_vals,TimeSeries(CycleOfInterest(i)*length(t_vals)+1:(CycleOfInterest(i)+1)*length(t_vals),5),'LineWidth',3,'Color',linecolors.I); %Plot I1
        semilogy(t_vals,TimeSeries(CycleOfInterest(i)*length(t_vals)+1:(CycleOfInterest(i)+1)*length(t_vals),7),'LineWidth',3,'Color',linecolors.L); %Plot L1
        semilogy(t_vals,TimeSeries(CycleOfInterest(i)*length(t_vals)+1:(CycleOfInterest(i)+1)*length(t_vals),9),'LineWidth',3,'Color',linecolors.V); %Plot V1
    
        semilogy(t_vals,TimeSeries(CycleOfInterest(i)*length(t_vals)+1:(CycleOfInterest(i)+1)*length(t_vals),4),'LineWidth',3,'Color',linecolors.E,'LineStyle','--'); %Plot E2  
        semilogy(t_vals,TimeSeries(CycleOfInterest(i)*length(t_vals)+1:(CycleOfInterest(i)+1)*length(t_vals),6),'LineWidth',3,'Color',linecolors.I,'LineStyle','--'); %Plot I2    
        semilogy(t_vals,TimeSeries(CycleOfInterest(i)*length(t_vals)+1:(CycleOfInterest(i)+1)*length(t_vals),8),'LineWidth',3,'Color',linecolors.L,'LineStyle','--'); %Plot L2    
        semilogy(t_vals,TimeSeries(CycleOfInterest(i)*length(t_vals)+1:(CycleOfInterest(i)+1)*length(t_vals),10),'LineWidth',3,'Color',linecolors.V,'LineStyle','--'); %Plot V2
        
        ylim([1e-2 1e10]);
        xlim([0 max(t_vals)]);
        xticks(0:4:24);
        set(gca,'YMinorTick','off','Box','off');
        %title('Total population');
        set(gca,'YTick',logspace(-2,10,5));
        set(gca,'FontSize',14,'FontWeight','bold','FontName','Times New Roman');
        xlabel(t,'Time (hr)','FontSize',14,"FontWeight",'bold','FontName','Times New Roman','Interpreter','latex');
        ylabel('Density (mL$^{-1}$)','FontSize',14,"FontWeight",'bold','FontName','Times New Roman');
        %ylabel(t,'Concentration (mL^{-1})','FontSize',14,"FontWeight",'bold');

   end
   

   %% Add plot labels
   nexttile(1);
   text(-30,1e10,'(A)','FontSize',14,'FontWeight','bold');
   nexttile(3);
   text(-7,1e10,'(B)','FontSize',14,'FontWeight','bold');
   nexttile(4);
   text(-7,1e10,'(C)','FontSize',14,'FontWeight','bold');
   nexttile(5);
   text(-7,1e10,'(D)','FontSize',14,'FontWeight','bold');
   nexttile(6);
   text(-7,1e10,'(E)','FontSize',14,'FontWeight','bold');

   plotlabel = ['B','C','D','E'];
   nexttile(1);
   for i = 1:4
       rectangle('Position',[CycleOfInterest(i)-1 1e-2 1 1e12],'FaceColor',[.8 .8 .8],'EdgeColor','none');
       text(CycleOfInterest(i),1e9,['$\rightarrow$ (' plotlabel(i) ')'],'FontSize',14,'FontWeight','bold');
   end
   
   %% Save Figure
filename = dir('..\Figures\FigureS4*');

if isempty(filename)
    filename = '..\Figures\FigureS4_v1.eps';
else
    filename = filename(end).name;
    version = extractBetween(filename,"_v",".");
    version = version{1};
    version = str2num(version);
    filename = ['..\Figures\' extractBefore(filename,num2str(version)) num2str(version+1) '.eps'];
end
filename1 = [filename(1:end-4) '.png'];

exportgraphics(fig,filename,"BackgroundColor",'none','ContentType','vector');
exportgraphics(fig,filename1,"BackgroundColor",'white');


   
