function [tideGauges, boundaryNodes] = tide_gauge_boundary_node_comparison(year, tideGaugePath, boundaryNodePath)
% The function reads Z0 data extracted from NCOM divided by into months for a set of nodes nearly latitudinally
% equivalent to a corresponding set of tide guages. 
%
% Assumed NOAA file name struct: 
% 	name_month_year.txt
%	neah_02_1999.txt
%
% Assumed NCOM derived Z0 file names:
%   Z0_month.th
%   Z0_02.th
%
% Input:
%   year			 - Year of data, assumes year is part of file path for tide data
%	tideGaugePath    - Path to directory holding tide guage data divided into monthly files
%	boundaryNodePath - Path to directory holding boudary node Z0 data from NCOM collected via ssh2d-format script
%
% Output:
% 	tideGauges    - Array of objects holding tide gauge station information
%	boundaryNodes - Array of object holding boundary nodes information 
%

% Ordered latitudinally from South to North - Made to match NOAA file names
stations = {'South Beach (Newport)', 'Garibaldi (Tillamook)', 'Tongue Point (Astoria)', 'Toke Point (Willapa)' ...
		    'Westport (Grays)', 'La Push', 'Neah Bay'};
fileStations = {'southbeach', 'garibaldi', 'tpoin', 'tokepoint', 'westport', 'lapush', 'neah'};
latitudes = (44.63333, 45.5533, 46.20000, 46.70666, 46.90333, 47.91333, 48.3666);


classdef station
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
		function obj = station(number, year, dataPath);
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

classdef boundNode
	% Boundary nodes class
	properties
		name		% Station name
		latitude 	% Latitude (decimal)
		number		% Station number. 1 being to most Southern, n being the most Northern
		monthlyZ0	% Monthly averaged Z0 - Initialized to NaNs
		year		% Year of this data
		dataPath	% Path to directory of data
		rawData		% Raw data from boundary nodes
	end

	methods
		function obj = boundaryNode(number, year, dataPath)
			obj.name 	  = station{number};
			obj.number	  = number;
			obj.lat  	  = latitudes(number);
			obj.year 	  = year;
			obj.dataPath  = dataPath;
			obj.monthlyZ0 = NaN(12);
		end
		function ret = calcZ0(obj, month)
			column = obj.number+2;
			systemCall = sprintf('awk ''{sum=sum+$%d} END {print sum/NR}'' Z0_%02d.th', column, month);
			[ret, obj.monthlyZ0(month)] = system(systemCall);
		end	
	end	
end

% For every station
%	Create object for that station
%   Create equivalent node boundary object 
% 	For every month 
%   	If there is tide gauge data
%			Load the month worth of NOAA tide gauge data
%			Calculate the monthly mean Z0 from the tide gauge data using UTide
%			Calculate the monthly mean Z0 from the nearest latitudinally equivalent node on the boundary
for statNumber = 1:size(stations,2)
	tideGauges(statNumber) = station(statNumber, year, tideGagePath); 
	boundaryNodes(statNumber) = boundaryNode(statNumber, year, boundaryNodePath); 
	for month = 1:12
		fileData = dir(tideGauges(statNumber).dataPath);
		if fileData.bytes ~= 0
			loadMonthData(tideGauges(statNumber, month);
			calcZ0(tideGauges(statNumber, month);
			calcZ0(boundaryNodes(statNumber, month);
		end
	end
end			 
	
