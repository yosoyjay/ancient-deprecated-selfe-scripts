classdef boundaryNode
	% Boundary nodes class
	properties
		name		% Station name
		lat	 		% Latitude (decimal)
		number		% Station number. 1 being to most Southern, n being the most Northern
		monthlyZ0	% Monthly averaged Z0 - Initialized to NaNs
		year		% Year of this data
		dataPath	% Path to directory of data
		rawData		% Raw data from boundary nodes
	end

	methods
		function obj = boundaryNode(number, name, fileName, latitude, year, dataPath);
			obj.name 	  = name;
			obj.number	  = number;
			obj.lat  	  = latitude;
			obj.year 	  = year;
			obj.dataPath  = dataPath;
			obj.monthlyZ0 = NaN(1,12);
		end
		function obj = calcZ0(obj, month)
			column = obj.number+2;
			systemCall = sprintf('awk ''{sum=sum+$%d} END {print sum/NR}'' %s/Z0_%02d.th', column, obj.dataPath, month);
			[ret, temp] = system(systemCall);
			obj.monthlyZ0(month) = str2num(temp);
		end	
	end	
end


