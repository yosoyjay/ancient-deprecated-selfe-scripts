year = 2000;
monthName = {'January', 'February', 'March', 'April', 'May', 'June', 'July', ...
			 'August', 'September', 'October', 'November', 'December'};

% Loop over elvation data and derive indices into data to divide data by month
% This is used to divide plots into months below. 
plotIndexing = NaN(12,2);
i = 1;
j = 2;
while j <= size(elev,1)
	if elev(i,5) == elev(j,5)
		j = j+1;
	else
		plotIndexing(elev(i,5), 1) = i;
		plotIndexing(elev(i,5), 2) = j-1;
		i = j;
		j = j+1;
	end
end

ssh_lats = [latitudes(1) latitudes(3) latitudes(4) latitudes(7)];
% Loop over months, if plotIndexing has values plot it 
% This assumes that getNCOMMonthlyMeanSSH has already been run and the NCOM values
% are held in meanSSH.
for i = 1:12
	if ~isnan(plotIndexing(i,1));
		range = plotIndexing(i,1):plotIndexing(i,2); 
		plot(elev(range,1), elev(range,2), 'color', 'r', 'marker', 'o'); hold on;
%		plot(elev(range,1), elev(range,4), 'color', 'b', 'marker', 'x'); hold on;
		plot(ssh_lats, meanSSH(:,i), 'color', 'b', 'marker', 'x'); hold on;
%		plot(elev(range,1), elev(range,2)-elev(range,4), 'color', 'k'); hold on;
		legend('Tide Gauge', 'Model');
		ylabel('Elevation m');
		xlabel('Latitude');
		title(monthName{i});
		figureName = sprintf('plots/%02d_%04d_ncom_tg.png', i, year); 
		saveas(gcf, figureName);				
		figureName = sprintf('plots/%02d_%04d_ncom_tg.fig', i, year); 
		saveas(gcf, figureName);				
		clf;
	end
end
