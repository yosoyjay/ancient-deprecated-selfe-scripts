function [tgStations] = tg_model_comparions(year, tidePath)
%  This script extracts data from tide gauges and boundary nodes.
%  Input:
%   year - int year
%   tidePath - Path the directory with tide gauge data

%  For every station 
%       Create object for station
%       Create object for boundary node
%       For every month 
%           If there is station data
%               Calculate elev for tide guage
%               Calculate elev for corresponding boundary node
load '/home/workspace/users/lopezj/data/tide_gauges/tg_stations.mat'
for station = 1:size(stations,2)
    tgStations(station) = tgStation(station, stations{station}, fileStations{station}, ... 
                                    latitudes(station), year, tidePath); 
	tic
    for month = 1:12
        filePath = sprintf('%s/%s_%d_%d.txt', tidePath, tgStations(station).fileName, month, year);
        fileData = dir(filePath);
        if (fileData.bytes ~= 0)
            tgStations(station) = tgStations(station).loadMonthData(month);
            tgStations(station) = tgStations(station).calcElev(month);
        end
    end
	toc
end
elev = [];

%index = 1;
%for month = 1:12
%   for station = 1:size(stations,2)
%       if ~isnan(tgStations(station).monthlyZ0(month))
%           elev(index,:) = [tgStations(station).lat, tgStations(station).monthlyElev(month), ...
%                            modelStations(station).lat, modelStations(station).monthlyElev(month), month];
%           index = index + 1;
%       end
%   end
%end

