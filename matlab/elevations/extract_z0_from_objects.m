index = 1;
for month = 1:12
	for station = 1:size(stations,2)
		if ~isnan(tgStations(station).monthlyZ0(month))
%			plot(tgStations(station).lat, tgStations(station).monthlyZ0(month), 'color', 'b'); hold on;
%			plot(boundaryNodes(station).lat, boundaryNodes(station).monthlyZ0(month), 'color', 'r'); hold on;
			z0(index,:) = [tgStations(station).lat, tgStations(station).monthlyZ0(month), ...
						 boundaryNodes(station).lat, boundaryNodes(station).monthlyZ0(month), month];
			index = index + 1;
		end
	end
end
