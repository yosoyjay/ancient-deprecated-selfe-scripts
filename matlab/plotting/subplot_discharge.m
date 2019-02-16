function [dischargePlot] = subplot_discharge(startTime, endTime)
% Creates a time series plot of the discharge at Beaver Army
%
% Input:
%	startTime - The starting time of the plot (datenum format)
%	endTime   - The ending time of the plot (datenum format)
%
% Output:
%  	dischargePlot - A handle to the time series figure (figure handle)
%
% lopezj - 01/25/2012
%

% Get the data
vars  = {'flux'};
table = {'external.discharge'};
startTime = datestr(startTime, 'mmmm dd, yyyy');
endTime   = datestr(endTime, 'mmmm dd, yyyy');
order = {'time'};
[discharge] = get_obs_data(vars, table, startTime, endTime, order);

% Plot the data
dischargePlot = figure('visible', 'off', 'color', 'w');
plot(discharge(:,1), discharge(:,2), ...
     'LineWidth', 1.5, ...
	 'Color', 'blue')

% Annotate the plot
set(gca,'FontSize', 14);
datetick('x');
set(gca,'xlim', [startTime endTime]); 
ylabel('Discharge m^3/s');
grid on;
ts = sprintf('Discharge at Beaver Army Terminal');
title(ts); 

% Bye
