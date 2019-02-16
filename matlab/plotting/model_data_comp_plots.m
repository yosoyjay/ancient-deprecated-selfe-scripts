function [] = model_data_comp_plots(model, obsInt, sa, variable, names) 
% [] = model_data_comp_plots(model, obsInt, sa, variable, names) 
% Creates a set of time serios plots comparing model and observation data.
%
% A unique plot created for every station and depth.
%
% Input:
%   model  - Model data returned from model_data_comp.m
%   obsInt - Observation data interpolated to model time step from model_data_comp.m
%   sa     - Skill assessment from model_data_comp.m
%	variable - The variable being plotted
%	names - Names for each run (cell array of strings)
%
% Output:
%   A set of plot saved as *.png
%
% lopezj - 01/22/2012
%

% varycolor
addpath '/home/workspace/users/lopezj/scripts/matlab/plotting';

% Map variables to ylabels
ylabels = containers.Map({'salt.63',      'temp.63',       'turb.63',       'elev.61'}, ...
                         {'Salinity psu', 'Temperature C', 'Turbidity NTU', 'Elevation m'});

% Plot set crap
%plotSet(1).color = 'red';
%plotSet(1).style = '--';


% Make a plots folder if it doesn't exists for the output 
plotsDir = sprintf('plots');
if (~exist(plotsDir,'file'))
	mkdir(plotsDir);
end

% Linewidth issues
if variable == 'elev.61'
	lineW = 1;
else
	lineW = 1.5;
end

% Loop over each station, depth, and model result to create plots
statNames = fieldnames(obsInt);
for station = 1:size(statNames,1)
    obsDepthNames = fieldnames(obsInt.(statNames{station}));
    for depth = 1:size(obsDepthNames,1)
        figH = figure('visible', 'off','color','w');
        % Plot observation data first
        plot(obsInt.(statNames{station}).(obsDepthNames{depth})(:,1), ... 
             obsInt.(statNames{station}).(obsDepthNames{depth})(:,2), ... 
             'LineWidth', lineW,     ...
             'Color', 'black');
        hold on;
        % Plot model data next
        runNames = fieldnames(model);
		plotSet = varycolor(length(runNames));
		modDepthNames = fieldnames(model.(runNames{1}).(statNames{station}));
        for run = 1:size(runNames,1) 
                plot(model.(runNames{run}).(statNames{station}).(modDepthNames{depth})(:,1), ...
                model.(runNames{run}).(statNames{station}).(modDepthNames{depth})(:,2), ...
                'LineWidth', lineW,  ...
                'Color', plotSet(run,:));
%                 'LineStyle', plotSet(run).style);
        end

        % Annotate plots
        set(gca,'FontSize', 14);
		xmin = model.(runNames{run}).(statNames{station}).(modDepthNames{1})(1,1);
		xmax = model.(runNames{run}).(statNames{station}).(modDepthNames{1})(end,1);
        datetick('x');
		set(gca,'xlim', [xmin xmax]); 
        ylabel(ylabels(variable));
		grid on;
        legStr = {'''Observation'''};
        for run = 1:size(runNames,1)
            modStr = sprintf(',''%s\tRMSE: %0.3f''', names{run}, ...
                             sa.(runNames{run}).(statNames{station}).(modDepthNames{depth})(6));
            legStr = strcat(legStr,modStr);
        end
        modStr = sprintf(',''Location'', ''NorthEastOutside''');
        legStr = strcat(legStr,modStr);
        eval(['l = legend(' legStr{1} ');']);
		set(l,'Interpreter','none');
		%screenSize = get(0,'ScreenSize')
		%set(l,'Position', [screenSize(3)/2 screenSize(4)*0.1 screenSize(3)/2 screenSize(4)*0.3]);
        dMeters = str2num(modDepthNames{depth}(5:end))/100;
        title(sprintf('%s %0.1f m Depth', statNames{station}, dMeters)); 

        % Save the plots
        fn = sprintf('%s/%s_%s_%0.1f_%s_%s.png', plotsDir, variable(1:4), statNames{station}, dMeters, datestr(xmin,'dd'), datestr(xmax,'dd'))
		iw = 1024;
		ih = 800;
		set(gcf,'PaperUnits','inches','PaperPosition',[0 0 iw/100.0 ih/100.0])
		print('-dpng', fn, '-r100');
		close(figH);
    end
end
