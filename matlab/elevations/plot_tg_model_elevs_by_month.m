year = 2000;
monthName = {'January', 'February', 'March', 'April', 'May', 'June', 'July', ...
			 'August', 'September', 'October', 'November', 'December'};

% Loop over elvation data and derive indices into data to divide data by month
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

% Loop over months, if plotIndexing has values plot it 
for i = 1:12
	if ~isnan(plotIndexing(i,1));
		range = plotIndexing(i,1):plotIndexing(i,2); 
		plot(elev(range,1), elev(range,2), 'color', 'r', 'marker', 'o'); hold on;
		plot(elev(range,1), elev(range,4), 'color', 'b', 'marker', 'x'); hold on;
		plot(ssh_lats, meanSSH(:,i), 'color', 'k', 'marker', '*'); hold on;
%		plot(elev(range,1), elev(range,2)-elev(range,4), 'color', 'k'); hold on;
		legend('Tide Gauge', 'DB22', 'NCOM');
		ylabel('Elevation m');
		xlabel('Latitude');
		title(monthName{i});
		figureName = sprintf('plots/%02d_%04d_month.png', i, year); 
		saveas(gcf, figureName);				
		figureName = sprintf('plots/%02d_%04d_moth.fig', i, year); 
		saveas(gcf, figureName);				
		clf;
	end
end
