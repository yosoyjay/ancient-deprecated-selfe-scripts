classdef tgStation
	% Tide gauge station class with no data hiding!
	properties
		name		% Station name
		fileName    % Name used as part of NOAA file
		number		% Station number. 1 being to most Southern, n being the most Northern
		lat	 		% Latitude (decimal)
		monthlyZ0	% Monthly averaged Z0
		monthlyElev % Monthly averaged elevation
		year		% Year of this data
		dataPath	% Path to directory of data
		rawData		% Raw data from NOAA tide gauge files, use as temp storage
	end

	methods
		function obj = tgStation(number, name, fileName, latitude, year, dataPath);
			obj.name 	  = name;
			obj.fileName  = fileName;
			obj.number    = number; 
			obj.lat  	  = latitude;
			obj.year 	  = year;
			obj.dataPath  = dataPath;
			obj.monthlyZ0 = NaN(1,12);
			obj.monthlyElev = NaN(1,12);
		end
		function obj = loadMonthData(obj, month)
			filePath = sprintf('%s/%s_%d_%04d.txt', obj.dataPath, obj.fileName, month, obj.year);
			fid = fopen(filePath, 'r');
			obj.rawData= textscan(fid, '%s %f', 'Delimiter', ',');  
			fclose(fid);
			% Convert time to datenum for use with UTide package
			obj.rawData{1,1} = datenum(obj.rawData{1,1});
		end
		function obj = calcZ0(obj, month)
			ret = ut_solv(obj.rawData{1,1}, obj.rawData{1,2}, [], obj.lat, 'auto'); 
			obj.monthlyZ0(month) = ret.mean;
		end				
		function obj = calcElev(obj, month)
			obj.monthlyElev(month) = mean(obj.rawData{1,2});
		end	
	end	
end

