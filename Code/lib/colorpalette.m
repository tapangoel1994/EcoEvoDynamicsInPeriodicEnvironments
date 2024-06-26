%% This script generates the marker and color palatte used for timeseries plots across this project

linecolors.R = .5*[1 1 1];
linecolors.S = '#0000ff';
linecolors.E = '#ffa600';
linecolors.I = '#a74ac7';
linecolors.L = '#00cbff';
linecolors.V = '#ff0000';

marker.lytic = 'd';
marker.temperate = 's';
marker.lysogenic = 'o';

linestyle.lytic = '-';
linestyle.temperate = '-.';
linestyle.lysogenic = ':';

% Set default interpreters to 'latex' for consistent font styling
% This script changes all interpreters from tex to latex. 

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));

for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end
% Optionally, set a consistent font style, e.g., 'Times New Roman'
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');



save('colorpalette.mat');