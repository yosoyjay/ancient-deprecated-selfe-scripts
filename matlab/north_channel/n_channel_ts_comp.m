function [model, obs, obsInt, sa] = n_channel_ts_comp(runDirs, startDay, endDay, file, stats)
% Script extracts model and observation data for the time given and returns the
% model and observation data as well as a comparison time series plot for each 
% station.
%
% Stations are hard coded based on my interests.
%
% Input:
%   runDirs  - List of run directories to compare multiple runs in a cell array
%            - '/path/to/run/runName/run' - runName will be used in struct if > 1 
%   startDay - Integer of the output file day. (1_salt.63, 8_hvel.64, etc.)
%   endDay   - Integer of the output file day.
%   file     - Variable of interest, must be *.63 ( 'salt.63', 'temp.63', etc.)
%	stats    - List of stations to compare model/observation data in a cell array
%
% Ouput:
%   model  - Model data for each station and time step
%   obs    - Observation data for each station and depth
%   obsInt - Observation data interpolated to model time
%   sa     - Skill assessment for every run, station, depth
%
% Struct: Matches that returned by get_model_station_data.m
%   obs.station.depth   - station corresponds to db name
%   obs.saturn2.dep_100 - Data at Saturn02 at depth of 1m
%   obs.saturn2.dep_600 - Data at Saturn02 at depth of 6m
%   obs.sandi.dep_490   - Data at Sandi at depth of 4.9m
%
%   model.baseline.saturn2.dep_100 - Model data of run baseline (among >1) at Saturn02 depth 1m
%       - If there is more than 1 run directory in runDirs and one of the dirs is 'baseline/run/'
%   model.run.saturn02.dep_100 - Model data of 1 run at Saturn02 depth 1m
%       - If there is only 1 run directory in runDirs
%
% lopezj - 01/17/2012
%
tic
% Paths and other constant stuff 
addpath '/home/workspace/users/lopezj/scripts/matlab/model_extraction/';
addpath '/home/workspace/users/lopezj/scripts/matlab/obs_extraction/';
addpath '/home/workspace/users/lopezj/scripts/matlab/skill_assessment/';

% Get model data using Joseph's Fortran code 
if length(runDirs) == 1
    model.run = get_model_station_data(runDirs{1}, startDay, endDay, stats, file, 0);
else
    for run = 1:length(runDirs)
        idx = findstr(runDirs{run},'/run')-1;
        names(run) = {runDirs{run}(1:idx)}
        model.(names{run}) = get_model_station_data(runDirs{run}, startDay, endDay, stats, file);
    end
end

% Calc model timestep and setup times for observation data
% endTime gets a timestep added to it so interpolation works correctly
if length(runDirs) == 1
    depths = fieldnames(model.run.(stats{1}));
    dt = (model.run.(stats{1}).(depths{1})(2,1) - model.run.(stats{1}).(depths{1})(1,1))*datenum(1)*86400;
else
    depths = fieldnames(model.(names{1}).(stats{1}));
    dt = (model.(names{1}).(stats{1}).(depths{1})(2,1) - model.(names{1}).(stats{1}).(depths{1})(1,1))*datenum(1)*86400;
end
startTime = get_start_date(runDirs{1})
endTime = addtodate(startTime, startDay-endDay+1, 'day')
endTime = addtodate(endTime, dt, 'second')
datestr(startTime)
datestr(endTime)

% Get observation data 
[obs, obsInt] = get_obs_station_data(startTime, endTime, stats, file, dt,1);

% Calculate skill assessment
if length(runDirs) == 1 
    for station = 1:length(stats)
        depths = fieldnames(obsInt.(stats{station}));
        for depth = 1:size(depths)
            temp = sa_combine(obsInt.(stats{station}).(depths{depth})(:,2), ...
                              model.run.(stats{station}).(depths{depth})(:,2));
            sa.run.(stats{station}).(depths{depth}) = temp;
        end
    end
else
    for run = 1:length(names)
        for station = 1:length(stats)
            depths = fieldnames(obsInt.(stats{station}));
            for depth = 1:size(depths)
                temp = sa_combine(obsInt.(stats{station}).(depths{depth})(:,2), ...
                                  model.(names{run}).(stats{station}).(depths{depth})(:,2));
                sa.(names{run}).(stats{station}).(depths{depth}) = temp;
            end
        end
    end
end
toc
