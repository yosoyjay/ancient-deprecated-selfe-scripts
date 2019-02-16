function [modelData] = get_model_station_data(runDir, startDay, endDay, stationList, variable, official)
% [modelData] = get_model_station_data(runDir, startDay, endDay, stationList, variable)
%
% Script extracts model data for each station listed in the station list which
% is a cell array.
%
% Input:
%   runDir   - The run directory
%   startDay - Integer of the output file day. (1_salt.63, 8_hvel.64, etc.)
%   endDay   - Integer of the output file day.
%	stationList - A cell array of stations to extract the data.
%	variable    - The model output variable to be extracted (e.g. salt.63, temp.63, etc..)
%	official - 1 = runDir to db2?, 0 = runDir to experimental run
%
% Output:
%	modelData 	- A structure that holds the extracted data
%
% Struct: Matches that returned by get_obs_station_data.m
%
% lopezj - 01/20/2012
%

% Paths and other constant stuff
addpath '/home/workspace/users/lopezj/scripts/matlab';
binPath     = '/home/workspace/users/lopezj/bin/selfe_post_processing';
stationPath = '/home/workspace/project/lopezj/data/stations';
fixedDepth  = 'read_output7b_group_z2';
freeSurface = 'read_output7b_group_zfs2';
eval(sprintf('load %s/station_map.mat', stationPath));
unitType = 1; % Everything in spcs

% Save cwd and change to run output directory
old_wd = pwd;
if official
	eval(sprintf('cd %s', runDir));
else
	eval(sprintf('cd %s/outputs', runDir));
end
%eval(sprintf('cd %s', runDir));

% Dates to datenum and interpolated to timestep of model
% Get dt and nSteps from '1_elev.61'
hElev = sz_readHeader('1_elev.61');
dt = hElev.dt;
nSteps = hElev.nSteps;

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

% Used to get elevations - one off hack.
if class(stationList) == 'char';
	% Create input file for extraction
	modelInPath = sprintf('model.in',runDir);
	modelIn = fopen(modelInPath,'w');
	text = sprintf('%d\n', 1);
	fwrite(modelIn, text);
	text = sprintf('%d\n', unitType);
	fwrite(modelIn, text);
	text = sprintf('%s\n', 'elev.61');
	fwrite(modelIn, text);
	text = sprintf('%d %d\n', startDay, endDay);
	fwrite(modelIn, text);
	fwrite(modelIn, '0');
	fclose(modelIn);

	% Weaksauce - code dies if the node is dry...
	eval(sprintf('!ln -fs %s/%s ./station.bp', stationPath, stationList)); 
	eval(sprintf('!%s/%s >> ../model_extraction.log 2>&1', binPath, fixedDepth));
	if official
		dataPath = sprintf('%s/fort.18',runDir);
	else
		dataPath = sprintf('%s/outputs/fort.18',runDir);
	end
	eval(sprintf('load %s', dataPath));

	depStr = sprintf('dep_0');
	modelData(:,1) = times';
	modelData(1:size(fort(:,2),1),2) = fort(:,2);
	return
end

% Loop over each station, get the data, and organize
for station = 1:length(stationList)
	% Get station information 
	statInfo = stationInfo(stationList{station});

	% Create input file for extraction
	modelInPath = sprintf('model.in',runDir);
	modelIn = fopen(modelInPath,'w');
	text = sprintf('%d\n', statInfo.stationType);
	fwrite(modelIn, text);
	text = sprintf('%d\n', unitType);
	fwrite(modelIn, text);
	text = sprintf('%s\n', variable);
	fwrite(modelIn, text);
	text = sprintf('%d %d\n', startDay, endDay);
	fwrite(modelIn, text);
	fwrite(modelIn, '0');
	fclose(modelIn);

	% Get the data
	% Depth from MSL 
	if statInfo.depthType == 1
		% Weaksauce - code dies if the node is dry...
		eval(sprintf('!ln -fs %s/%s ./station.sta', stationPath, statInfo.stationFile)); 
		eval(sprintf('!%s/%s >> ../model_extraction.log 2>&1', binPath, fixedDepth));
		if official
			dataPath = sprintf('%s/fort.18',runDir);
		else
			dataPath = sprintf('./fort.18');
		end
		eval(sprintf('load %s', dataPath));
	    %modelData.(stationList{station}) = fort;
		
		% Organize data
		for depth = 1:size(statInfo.depths,2)
			depStr = sprintf('dep_%i', int16(statInfo.depths(depth)*100));
			modelData.(stationList{station}).(depStr)(:,1) = times';
			modelData.(stationList{station}).(depStr)(1:size(fort(:,depth+1)),2) = fort(:,depth+1);
		end
	% Depth from free surface
	elseif statInfo.depthType == 2
		% Weaksauce - code dies if the node is dry...
		eval(sprintf('!ln -fs %s/%s ./station.sta', stationPath, statInfo.stationFile)); 
		eval(sprintf('!%s/%s >> ../model_extraction.log 2>&1', binPath, freeSurface));
		if official
			dataPath = sprintf('%s/fort.18',runDir);
		else
			dataPath = sprintf('./fort.18');
		end
		eval(sprintf('load %s', dataPath));
	    %modelData.(stationList{station}) = fort;
		
		% Organize data
		for depth = 1:size(statInfo.depths,2)
			depStr = sprintf('dep_%i', int16(statInfo.depths(depth)*100));
			modelData.(stationList{station}).(depStr)(:,1) = times';
			modelData.(stationList{station}).(depStr)(1:size(fort(:,depth+1)),2) = fort(:,depth+1);
		end
	% Used to get elevations 
	else
		% Weaksauce - code dies if the node is dry...
		eval(sprintf('!ln -fs %s/%s ./station.bp', stationPath, statInfo.stationFile)); 
		eval(sprintf('!%s/%s >> ../model_extraction.log 2>&1', binPath, fixedDepth));
		if official
			dataPath = sprintf('%s/fort.18',runDir);
		else
			dataPath = sprintf('./fort.18');
		end
%		dataPath = sprintf('%s/fort.18',runDir);
		eval(sprintf('load %s', dataPath));
		
		% Organize data
		for depth = 1:size(statInfo.depths,2)
			depStr = sprintf('dep_%i', int16(statInfo.depths(depth)*100));
			modelData.(stationList{station}).(depStr)(:,1) = times';
			modelData.(stationList{station}).(depStr)(1:size(fort(:,depth+1),1),2) = fort(:,depth+1);
		end

	end
end
	
% Back to original directory
eval(sprintf('cd %s', old_wd));

