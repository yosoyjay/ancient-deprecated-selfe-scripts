function [] = model_comp_plots(model, variable, names) 
% Creates a set of time series plots comparing model and observation data.
%
% A unique plot created for every station and depth.
%
% Input:
%   model  - Model data returned from model_data_comp.m
%	variable - The variable being plotted
%	names - Names for each run (cell array of strings)
%
% Output:
%   A set of plot saved as *.png
%
% lopezj - 02/11/2012
%

% Map variables to ylabels
ylabels = containers.Map({'salt.63',      'temp.63',       'turb.63',       'elev.61'}, ...
                         {'Salinity psu', 'Temperature C', 'Turbidity NTU', 'Elevation m'});

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

% Linewidth issues
if variable == 'elev.61'
	lineW = 1;
else
	lineW = 1.5;
end

% Loop over each station, depth, and model result to create plots
runNames = fieldnames(model);
statNames = fieldnames(model.(runNames{1}));
for station = 1:size(statNames,1)
    modDepthNames = fieldnames(model.(runNames{1}).(statNames{station}));
    for depth = 1:size(modDepthNames,1)
        figH = figure('visible', 'off','color','w');
		hold on;
        % Plot model data next
        for run = 1:size(runNames,1) 
            plot(model.(runNames{run}).(statNames{station}).(modDepthNames{depth})(:,1), ...
                 model.(runNames{run}).(statNames{station}).(modDepthNames{depth})(:,2), ...
                 'LineWidth', lineW,  ...
                 'Color', plotSet(run).color, ...
                 'LineStyle', plotSet(run).style);
        end

        % Annotate plots
        set(gca,'FontSize', 14);
		xmin = model.(runNames{run}).(statNames{station}).(modDepthNames{1})(1,1);
		xmax = model.(runNames{run}).(statNames{station}).(modDepthNames{1})(end,1);
        datetick('x');
		set(gca,'xlim', [xmin xmax]); 
        ylabel(ylabels(variable));
		grid on;
%        legStr = {''};
%        for run = 1:size(runNames,1)
%            modStr = sprintf(',''%s''', names{run});
%            legStr = strcat(legStr,modStr);
%        end
%        modStr = sprintf(',''Location'', ''NorthEastOutside''');
%        legStr = strcat(legStr,modStr);
%        eval(['l = legend(' legStr{1} ');']);
%		set(l,'Interpreter','none');
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
