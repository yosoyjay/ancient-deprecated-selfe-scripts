function [tgStations, modelStations_14, modelStations_16, modelStations_22, elev] = tg_model_comparions(year, tidePath)
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
    modelStations_14(station) = modelStation(station, stations{station},  ... 
                                          latitudes(station), longitudes(station), 14, year);
    modelStations_16(station) = modelStation(station, stations{station},  ... 
                                          latitudes(station), longitudes(station), 16, year);
    modelStations_22(station) = modelStation(station, stations{station},  ... 
                                          latitudes(station), longitudes(station), 22, year);
	tic
    for month = 1:12
        filePath = sprintf('%s/%s_%d_%d.txt', tidePath, tgStations(station).fileName, month, year);
        fileData = dir(filePath);
        if (fileData.bytes ~= 0)
            tgStations(station) = tgStations(station).loadMonthData(month);
            tgStations(station) = tgStations(station).calcElev(month);
        end
        modelStations_14(station) = modelStations_14(station).calcElev(month);
        modelStations_16(station) = modelStations_16(station).calcElev(month);
        modelStations_22(station) = modelStations_22(station).calcElev(month);
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

