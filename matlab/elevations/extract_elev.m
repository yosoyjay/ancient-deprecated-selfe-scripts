index = 1;
for month = 1:12
	for station = 1:size(stations,2)
		if ~isnan(tgStations(station).monthlyElev(month))
			elev(index,:) = [tgStations(station).lat, tgStations(station).monthlyElev(month), ...
						 	 modelStations(station).lat, modelStations(station).monthlyElev(month), month];
			index = index + 1;
		end
	end
end

