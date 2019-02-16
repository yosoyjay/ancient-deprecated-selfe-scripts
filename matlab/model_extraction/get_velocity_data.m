function [velData, szDepths] = get_velocity_data(runDir, startDay, endDay, location)
% [velData] = get_velocity_data(runDir, startDay, endDay, location)
%
% This script extracts the velocity data at given location and returns it.
%
% Input:
%	runDir   - The run directory
%   startDay - Integer of the output file day. (1_salt.63, 8_hvel.64, etc.)
%   endDay   - Integer of the output file day.
%	location - X and Y of the location in vector [x y];
%
% Output:
%	velData	 - Velocity data at all nodes in the vertical (day,nSteps, nNodes, [uVel vVel wVel])
%
% lopezj - 03/21/2012
%

% Paths and other constant stuff
addpath '/usr/local/cmop/matlab/cmop/m-elio';
addpath '/home/workspace/users/lopezj/scripts/matlab/';

% Save cwd and change to run output directory
old_wd = pwd;
if official
	eval(sprintf('cd %s', runDir));
else
	eval(sprintf('cd %s/outputs', runDir));
end

% Read elev file to get information about the run
% Dates to datenum and interpolated to timestep of model
% Get dt, nSteps, and nVertLevels from '1_salt.63',
hSalt = sz_readHeader('1_salt.63');
dt = hSalt.dt;
nSteps = hSalt.nSteps;
nLevels = hSalt.vgrid.nLevels; 
nNodes = hSalt.hgrid.np;

% Determine time that correlates with files given based on data in bctides.in
% Assume that each output file is equal to 1 day
startDate = get_start_date('../');
%startDate = datenum('August 13, 2010 00:00:00');
startTime = addtodate(startDate, startDay-1, 'day');
startTime = addtodate(startTime, dt, 'second'); 
endTime   = addtodate(startDate, endDay, 'day');
times = startTime:1/nSteps:endTime;
datestr(startTime);
datestr(endTime);
size(times);

% Allocate space for velocity data (days, timeStep, node, level, [uVel vVel wVel])
velData = nans(endDay-startDay+1,nSteps,nNodes,nLevels,3);

% Loop over each day and extract the data
for day = startDay:endDay
	% Read headers of velocity files
	hvelPath = sprintf('%d_hvel.64', day);
	hHvel = sz_readHeader(hvelPath);
	vertPath = sprintf('%d_vert.63', day);
	hVert = sz_readHeader(vertPath);

	% Get data for all time steps
	[hvelData ts] = sz_readTimeStep(hHvel, 1:nSteps);
	[vertData ts] = sz_readTimeStep(hVert, 1:nSteps);

	% Map sz levels to depths - Use map_sz2hts_mat to do all timesteps at once
	dHvel = map_sz2hts_mat(hHvel, hvelData);
	dVert = map_sz2hts_mat(hVert, vertData);

	% Combine this data into one array and shove in (days, :, :, :, [uVel vVel wVel])	

