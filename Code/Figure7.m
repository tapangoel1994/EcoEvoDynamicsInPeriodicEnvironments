%%Code takes a range of one-host one-virus system with different (q,gamma) values and
%%plays each
%%strategy over 100 growth cycles each with randomly drawn passage
%%probabilities. This is repeated a 100 times to see in how many cases the strategy persists.
%%other parameters fixed. Uses an RSEILV model with no MOI dependence.


%% Date Created: 6/26/2024
%% Author: Tapan Goel

%% Life history parameters (units of hours, micrograms and mL). 
clear all;
close all;
rng(0);

addpath('utils\');
addpath('lib\');
    colorpalette;
    fixedparameters;


%% Additional simulation parameters
params.T = 24; % hours
params.t_vals = transpose(0:params.dt:params.T); % time
MaxCycles = 100;
NumSims = 1;

%% filter parameters
p_S = 0;
p_E = 0;
p_I = 0;
p_L = .1*ones(MaxCycles,NumSims);%*(rand(MaxCycles,NumSims));
p_V = 10.^(-5+4*(rand(MaxCycles,NumSims)));


%% initial conditions
R0 = 1e2; %initial resource amount in ug/mL ( 500 mL flask)
S0 = 1e7; %Initial concentration of susceptibles in flask (per mL)
V01= 1e4; %initial concentration of virus in flask (per mL)
V02 = 0;
x0 = [R0 S0 zeros(1,6) V01 V02];

QGamma = [0 0;.5 .083;1 0]; %% each row is the (q,gamma) for a particular strategy
StrategyLabels = {'Obligately lytic','Temperate','Obligately lysogenic'};

h = figure('Position',[100 0 600 1200]);
t = tiledlayout(4,1,"TileSpacing",'compact','Padding','compact');

%% Plot the viral fraction
nexttile(1);
semilogy(.99:MaxCycles-.01,p_V,'-*','Color',.25*[1 1 1],'LineWidth',2);
ylim([10^-5.1 1e-1]);
xlim([0 MaxCycles+5]);
set(gca,'YMinorTick','off','Box','off','XTick',0:20:MaxCycles,'XTickLabel',[],'YTick',logspace(-5,-1,3));
set(gca,'FontSize',14,'FontWeight','bold');
ylabel('$q_V$','FontSize',16,'FontWeight','bold');
 
%% Calculate and plot viral trajectories
for index = 1:3

    y0 = x0;
    params.q = [QGamma(index,1) 0];
    params.gamma = [QGamma(index,2) 0];

    timeseries = [];

    for iter = 1:MaxCycles

        [t_vals,y] = ode113(@ODE_RSEILV_2Species,params.t_vals,y0,options,params);
        timeseries = [timeseries;y(end,:)];
        TransferMatrix = diag([0 p_S p_E p_E p_I p_I p_L(iter,1) p_L(iter,1) p_V(iter,1) p_V(iter,1)]);
        y0 = [R0 S0 zeros(1,8)] + y(end,:)*TransferMatrix; 
    end

    nexttile(index+1);

    %semilogy(.99:MaxCycles,timeseries(:,1),'LineWidth',2,'Color',linecolors.R);
    
    semilogy(.99:MaxCycles,timeseries(:,2),'LineWidth',2,'Color',linecolors.S);
    hold on;
    semilogy(.99:MaxCycles,timeseries(:,3),'LineWidth',2,'Color',linecolors.E);
    semilogy(.99:MaxCycles,timeseries(:,5),'LineWidth',2,'Color',linecolors.I);
    semilogy(.99:MaxCycles,timeseries(:,7),'LineWidth',2,'Color',linecolors.L);
    semilogy(.99:MaxCycles,timeseries(:,9),'LineWidth',2,'Color',linecolors.V);

    %semilogy(.99:5:MaxCycles,timeseries(1:5:end,1),'LineWidth',1,'Color',linecolors.R);
    semilogy(.99:5:MaxCycles,timeseries(1:5:end,2),'LineWidth',1,'Color',linecolors.S,'Marker',marker.(Strategies{index}),'LineStyle','none');
    semilogy(.99:5:MaxCycles,timeseries(1:5:end,3),'LineWidth',1,'Color',linecolors.E,'Marker',marker.(Strategies{index}),'LineStyle','none');
    semilogy(.99:5:MaxCycles,timeseries(1:5:end,5),'LineWidth',1,'Color',linecolors.I,'Marker',marker.(Strategies{index}),'LineStyle','none');
    semilogy(.99:5:MaxCycles,timeseries(1:5:end,7),'LineWidth',1,'Color',linecolors.L,'Marker',marker.(Strategies{index}),'LineStyle','none');
    semilogy(.99:5:MaxCycles,timeseries(1:5:end,9),'LineWidth',1,'Color',linecolors.V,'Marker',marker.(Strategies{index}),'LineStyle','none');
    
    ylim([1e-2 1e10]);
    xlim([0 MaxCycles+5]);
    set(gca,'YMinorTick','off','Box','off','XTick',0:20:MaxCycles,'XTickLabel',[],'YTick',logspace(-2,10,5));
    text(50,1e9,StrategyLabels{index},'FontSize',14,'FontWeight','bold');
    if index == 3
        set(gca,'XTickLabel',0:20:MaxCycles);
    end
    set(gca,'FontSize',14,'FontWeight','bold');
    ylabel('Density (mL$^{-1}$)','FontSize',16,'FontWeight','bold');
    if index == 3
        xlabel('Number of cycles','FontSize',16,'FontWeight','bold')
    end
     

    if index == 1
        nexttile(2);
        lgd = legend('S','E','I','L','V');
        lgd.Position = [0.8165    0.5306    0.1590    0.1062];
        legend('boxoff');
    end

end

%% Add panel labels
nexttile(1);
text(-16,10^-1,'(A)','FontSize',16,'FontWeight','bold');
nexttile(2);
text(-16,10^10,'(B)','FontSize',16,'FontWeight','bold');
nexttile(3);
text(-16,10^10,'(C)','FontSize',16,'FontWeight','bold');
nexttile(4);
text(-16,10^10,'(D)','FontSize',16,'FontWeight','bold');
%% Save Figure
filename = dir('RandomSelections*');
filename = filename(end).name;
if isempty(filename)
    filename = 'SingleCycle_v1.eps';
else
    version = extractBetween(filename,"_v",".");
    version = version{1};
    version = str2num(version);
    filename = [extractBefore(filename,num2str(version)) num2str(version+1) '.eps'];
end
saveas(h,filename,'epsc');


