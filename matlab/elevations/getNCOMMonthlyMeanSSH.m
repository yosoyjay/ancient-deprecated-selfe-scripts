function [monthlyMeanSSH] = getNCOMMonthlyMeanSSH(year)
%  This script gets the monthly mean SSH for 4 stations:
%  South Beach, Tounge Point, Toke Point, and Neah Bay 
%
%  The indices for the locations are hardcoded and were not derived
%  programatically. Calls getNCOMSSH.m to get data from NCOM netCDF files. 
%
%  Input: 
%  	year - Made to only work with 2000 data YMMV
%
%  Output:
%	monthlyMeanSSH  - SSH for the 4 stations mentioned above averaged per month
%
% lopezj - 11/06/2011

% Hardcoded because that's what I'm looking at, plus I don't even use it!
% Take that memory!
year = 2000;

% Prep for looping and data formatting 
month = 1;
monthlyMeanSSH = NaN(4, 12);
monthly = NaN(4, 31);
date = sprintf('12-31-%04d', year-1);
date = datenum(date);

% Loop over each day of the year and get NCOM data for each station
% Separate into months and average 
% Results for each month go into monthlyMeanSSH by station 
for i = 1:367
	date = addtodate(date,1,'day');
	dv = datevec(date);
% If date has advanced to new month, average old and put in output array 
	if dv(2) ~= month
		monthlyMeanSSH(1,month) = mean(monthly(1,:));
		monthlyMeanSSH(2,month) = mean(monthly(2,:));
		monthlyMeanSSH(3,month) = mean(monthly(3,:));
		monthlyMeanSSH(4,month) = mean(monthly(4,:));
		month = month+1;
	end
% Get NCOM data and separate for data near tide gauges
	[tempSSH, lat, long] = getNCOMSSH(dv);
	monthly(1,dv(3)) = tempSSH(70, 295);	% South Beach
	monthly(2,dv(3)) = tempSSH(70, 304);	% Tongue Point
	monthly(3,dv(3)) = tempSSH(68, 320);	% Toke Point
	monthly(4,dv(3)) = tempSSH(63, 341);    % Neah Bay
end	
		
