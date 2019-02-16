function [ssh, lat, long] = getNCOMSSH(date)
% This script gets NCOM SSH for a given date.
% The best way I found to plot it is to use contourf(long,lat,ssh)
%
% I used datatip to get indices of specific nodes.
% 
% Input: 
%	date - datevec format
% Ouput:
%	ssh  - SSH data
%	lat  - Corresponding latitudes
%	long - Corresponding longitudes
%
% lopezj - 11/03/2011

% Base directory where all the NCOM data is
baseDir = '/home/workspace/ccalmr6/nrldata/';

% Deal with the dates and parse info for file formats and clarity
year  = date(1);
month = date(2);
day   = date(3);

% Get ssh data
sshFile = sprintf('%s%04d/ssh/ssh.glb8_2f_%04d%02d%02d00.nc', baseDir, year, year, month, day);
ncid    = netcdf.open(sshFile, 'NC_NOWRITE');
varid   = netcdf.inqVarID(ncid, 'X_Index');
X = netcdf.getVar(ncid, varid);
startX = X(1);
countX = size(X,1);
varid  = netcdf.inqVarID(ncid, 'Y_Index');
Y = netcdf.getVar(ncid, varid);
startY = Y(1);
countY = size(Y,1);
varid  = netcdf.inqVarID(ncid, 'Surface_Elevation');
ssh = netcdf.getVar(ncid, varid, 'double');
netcdf.close(ncid);

% Get lat/long associated with netcdf data
latFile = sprintf('%s%04d/model_lat.nc', baseDir, year);
ncid    = netcdf.open(latFile, 'NC_NOWRITE');
varid   = netcdf.inqVarID(ncid, 'Lat');
lat     = netcdf.getVar(ncid, varid, [startX startY], [countX countY]);
netcdf.close(ncid);

% Longitude is measured E of GWM, so it is adjust to reflect being W of GWM
longFile = sprintf('%s%04d/model_lon.nc', baseDir, year);
ncid     = netcdf.open(longFile, 'NC_NOWRITE');
varid    = netcdf.inqVarID(ncid, 'Long');
long     = netcdf.getVar(ncid, varid, [startX startY], [countX countY]);
long     = long-360;
netcdf.close(ncid);
