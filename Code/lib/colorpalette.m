%% This script generates the color scheme used for creating graphics accross this project.
%% Before running any code to generate graphs etc, run this script. 


%% Colors for R, S, E, I, L and V in line plots
linecolors.R = .5*[1 1 1];
linecolors.S = '#0000ff';
linecolors.E = '#ffa600';
linecolors.I = '#c568e5';%'#e386ff';%'#a74ac7';
linecolors.L = '#00cbff';
linecolors.V = '#ff0000';

%% Markers  and line styles for lytic, temperate and lysogenic strategies
marker.lytic = 'd';
marker.temperate = 's';
marker.lysogenic = 'o';

linestyle.lytic = '-';
linestyle.temperate = '-.';
linestyle.lysogenic = ':';
Strategies = {'lytic','temperate','lysogenic'};

%% Set default interpreters to 'latex' for consistent font styling
%% This script changes all interpreters from tex to latex. 

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
index_fontname = find(contains(list_factory,'FontName'));


for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

for i = 1:length(index_fontname)
    default_name = strrep(list_factory{index_fontname(i)},'factory','default');
    set(groot, default_name,'Times New Roman');
end

%% Optionally, set a consistent font style, e.g., 'Times New Roman'
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot,'defaultLineMarkerSize',8);

%% Colormap for PIP plots
PIPColorMap = [.4 .4 .4;... %invasion invalid because no resident
               1 1 1;... %invasion failed
               .8 .8 .8;...%resident =  mutant
               0 0 0]; % invasion success

%save('colorpalette.mat');