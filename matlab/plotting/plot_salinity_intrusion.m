function [silData] = plot_salinity_intrusion(silPaths, names, psu, startDay)
% function [silData] = plot_salinity_intrusion(silPaths, names, psu, startDay)
%
% Creates a time series plot comparing salinity intrusion length between a set of runs
%
% Input:
%   silPaths - A set of paths to extracted salinity lengths in a cell array
%   names    - Names of the runs used for plotting purposes
%   psu      - The psu value that defines the salinity intrusion length
%   startDay - Day to start the plot. (Integer day of run 1_.. 2_... etc...)
% Output:
%   intData(:,2) - Intrusion length data (time, intrusion length)
%
% lopezj - 1/24/2012
%

addpath '/home/workspace/users/lopezj/scripts/matlab/model_extraction/';

% Load the data
for run = 1:length(silPaths)
	% Get model run data
	runPath = sprintf('%s/../run', silPaths{run});
	runInfo = get_run_info(runPath);

	% Create for time for the run
	startTime = datenum(runInfo.startDate);
	endTime   = addtodate(startTime, runInfo.nDays, 'day');
	startTime = addtodate(startTime, runInfo.dt, 'second');
	startTime = addtodate(startTime, startDay-1, 'day');
	times     = startTime:1/runInfo.nSteps:endTime;

	% Load data from file and shove in struct
    file = sprintf('intrusion_length%i', psu);
    eval(['load ' silPaths{run} '/' file '.dat']);
    data = eval(file);
    silData.(names{run})(:,1) = times';
    silData.(names{run})(:,2) = data(:,4); 
end

% Plot set crap
plotSet(1).color = 'red';
plotSet(1).style = '--';
plotSet(2).color = 'blue';
plotSet(2).style = '-.';
plotSet(3).color = 'green';
plotSet(3).style = '--';
plotSet(4).color = 'cyan';
plotSet(4).style = '-.';
plotSet(5).color = 'black';
plotSet(5).style = '--';

% Make a plots folder if it doesn't exists for the output 
plotsDir = sprintf('plots');
if (~exist(plotsDir,'file'))
    mkdir(plotsDir);
end

% Create plot
figH = figure('visible', 'off','color','w');

% Loop over all run data
for run = 1:length(silPaths)
    % Plot first salinity intrusion 
    plot(silData.(names{run})(:,1), silData.(names{run})(:,2), ... 
         'LineWidth', 1.5,  ...
         'Color', plotSet(run).color, ...
         'LineStyle', plotSet(run).style);
    hold on;
end

% Annotate plots
set(gca,'FontSize', 14);
datetick('x');
set(gca,'xlim', [startTime endTime]); 
ylabel('Salinity intrusion [m]');
grid on;

% Legend crap
legStr = {};
for run = 1:length(names)
    tmp = sprintf('''%s'',', names{run});
    legStr = strcat(legStr,tmp);
end
legStr = legStr{1}(1:end-2);
legStr = strcat(legStr,''',''Location'',''NortheastOutside''');
eval(['l=legend(' legStr ');']);
set(l,'Interpreter', 'none');
ts = sprintf('Comparison of Salinity Intrusion (%i psu)', psu);
title(ts); 

% Save the plots
fn = sprintf('%s/sal_intrusion_%i.png', plotsDir, psu);
iw = 1024;
ih = 800;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 iw/100.0 ih/100.0])
print('-dpng', fn, '-r100');
close(figH);

