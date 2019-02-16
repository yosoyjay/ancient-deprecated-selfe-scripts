function [tgStations, boundaryNodes] = tg_bn_comparions(year, tidePath, boundNodePath)
%  This script extracts data from tide gauges and boundary nodes.

%  For every station 
%		Create object for station
%		Create object for boundary node
%		For every month 
%			If there is station data
%				Calculate Z0 for tide guage
%				Calculate Z0 for corresponding boundary node
load 'tg_stations.mat'
for station = 1:size(stations,2)
	tgStations(station) = tgStation(station, stations{station}, fileStations{station}, ... 
									latitudes(station), year, tidePath); 
	boundaryNodes(station) = boundaryNode(station, stations{station}, fileStations{station}, ... 
										  latitudes(station), year, boundNodePath);
	for month = 1:12
		filePath = sprintf('%s/%s_%d_%d.txt', tidePath, tgStations(station).fileName, month, year);
		fileData = dir(filePath);
		if (fileData.bytes ~= 0)
			tgStations(station) = tgStations(station).loadMonthData(month);
			tgStations(station) = tgStations(station).calcZ0(month);
			boundaryNodes(station) = boundaryNodes(station).calcZ0(month);
		end
	end
end

%for station = 1:size(stations,2)
%	for month = 1:12
%		if ~isnan(tgStations(station).monthlyZ0(month))
%			plot(tgStations(station).lat, tgStations(station).monthlyZ0(month)); hold on;
%			plot(boundaryNodes(station).lat, boundaryNodes(station).monthlyZ0(month)); hold on;
%		end
%	end
%end


