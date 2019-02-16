classdef tgStation
	% Tide gauge station class with no data hiding!
	properties
		name		% Station name
		fileName    % Name used as part of NOAA file
		number		% Station number. 1 being to most Southern, n being the most Northern
		latitude 	% Latitude (decimal)
		monthlyZ0	% Monthly averaged Z0
		year		% Year of this data
		dataPath	% Path to directory of data
		rawData		% Raw data from NOAA tide gauge files, use as temp storage
	end

	methods
		function obj = tgStation(number, year, dataPath);
			obj.name 	  = stations{number};
			obj.fileName  = fileStations{number};
			obj.number    = number; 
			obj.lat  	  = latitudes(number);
			obj.year 	  = year;
			obj.dataPath  = dataPath;
			obj.monthlyZ0 = NaN(12);
		end
		function loadMonthData(obj, month)
			filePath = sprintf('%s_%02d_%04d.txt', obj.name, month, obj.year);
			fid = fopen(filePath, 'r'); 
			obj.rawData = textscan(fid, '%s %f', 'Delimiter', ',');  
			fclose(fid);
			% Convert time to datenum for use with UTide package
			obj.rawData{1,1} = datenum(obj.rawData{1,1});
		end
		function calcZ0(obj, month)
			ret = ut_solv(obj.rawData{1,1}, obj.rawData{1,2}, [], obj.lat, 'auto'); 
			obj.monthlyZ0(month) = ret.mean;
		end				
	end	
end

