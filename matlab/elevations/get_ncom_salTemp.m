function [salt, temp, depth, lat, long] = get_ncom_salTemp(date)
% This script gets NCOM salt and temp for a given date.
% The best way I found to plot it is to use contourf(long,lat,ssh)
% 
% Input: 
%	date - datevec format
% Ouput:
%   salt[][][] - salt at every level (40 levels)
%   temp[][][] - temp at every level  
%   depth[] - Array of depths at the 40 levels for 
%	lat  - Corresponding latitudes
%	long - Corresponding longitudes
%
% lopezj - 01/04/2011
%

% Base directory where all the NCOM data is
baseDir = '/home/workspace/ccalmr6/nrldata/';

% Assume input date is datevec so no casting
year  = date(1);
month = date(2);
day   = date(3);

% Get salt and temp data 
% X & Y are the indices that correspond to each value. These are then used
% to get the corresponding lat and long values from different files.
saltFile = sprintf('%s%04d/s3d/s3d.glb8_2f_%04d%02d%02d00.nc', ...
                   baseDir, year, year, month, day);

tempFile = sprintf('%s%04d/t3d/t3d.glb8_2f_%04d%02d%02d00.nc', ...
                   baseDir, year, year, month, day);

% Salt Salt Salt Salt Salt Salt Salt FACE
% X index 
ncid    = netcdf.open(saltFile, 'NC_NOWRITE');
varid   = netcdf.inqVarID(ncid, 'X_Index');
X = netcdf.getVar(ncid, varid);
startX = X(1);
countX = size(X,1);
% Y index 
varid  = netcdf.inqVarID(ncid, 'Y_Index');
Y = netcdf.getVar(ncid, varid);
startY = Y(1);
countY = size(Y,1);
% Salt data 
varid  = netcdf.inqVarID(ncid, 'Salinity');
salt = netcdf.getVar(ncid, varid, 'double');
netcdf.close(ncid);

% Temp Temp Temp TEmp TEmp TEMp TEMP
% Temp data
ncid    = netcdf.open(tempFile, 'NC_NOWRITE');
varid  = netcdf.inqVarID(ncid, 'Temperature');
temp = netcdf.getVar(ncid, varid, 'double');
netcdf.close(ncid);

% Get lat/long associated with netcdf data
latFile = sprintf('%s/model_lat.nc', baseDir);
ncid    = netcdf.open(latFile, 'NC_NOWRITE');
varid   = netcdf.inqVarID(ncid, 'Lat');
lat     = netcdf.getVar(ncid, varid, [startX startY], [countX countY]);
netcdf.close(ncid);

% Longitude is measured E of GWM, so it is adjust to reflect being W of GWM
longFile = sprintf('%s/model_lon.nc', baseDir);
ncid     = netcdf.open(longFile, 'NC_NOWRITE');
varid    = netcdf.inqVarID(ncid, 'Long');
long     = netcdf.getVar(ncid, varid, [startX startY], [countX countY]);
long     = long-360;
netcdf.close(ncid);

% Get depths for each level - Z is hard coded because there is always 40 levels
depthFile = sprintf('%s/model_zm.nc', baseDir);
ncid     = netcdf.open(depthFile, 'NC_NOWRITE');
varid    = netcdf.inqVarID(ncid, 'zm');
depth    = netcdf.getVar(ncid, varid, [startX startY 1], [countX countY 39]);
netcdf.close(ncid);

