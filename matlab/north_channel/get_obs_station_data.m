function [obsData, obsInterp] = get_obs_station_data(startDay, endDay, stationList, variable, modelTs)
% Extracts and returns observation data that matches those list stationList cell array.
%
% Input:
%   startDay    - First day of data in datenum format
%   endDay      - End day of data in datenum format
%   stationList - Cell array of stations to get data
%   variable    - The variable to extract (e.g. 'salt.63', 'temp.63', etc...)
%   modelTs     - Model time step, used to interpolated observation data to model 
%                 time step for SA
%
% Output:
%   obsData   - A structure that holds the observation data
%   obsInterp - A structure that holds the observation data interpolated to 
%               model time steps
%
% Struct: Matches that returned by get_model_station_data.m
%   obsData.station.depth   - station corresponds to db name
%   obsData.saturn2.dep_100 - Data at Saturn02 at depth of 1m
%   obsData.saturn2.dep_600 - Data at Saturn02 at depth of 6m
%   obsData.sandi.dep_490   - Data at Sandi at depth of 4.9m
%
% lopezj - 01/22/2012
%

% Paths and other constant stuff
addpath '/home/workspace/users/lopezj/scripts/matlab/';

% Map variables to database term
varMap = containers.Map({'temp.63',     'salt.63',  'turb.63'},         ...
                        {'temperature', 'salinity', 'turbidity'});

% Connect to database
dbConn = connect_to_db();

% Seperate time and data in two different queries to minimize messing with Matlab
% formatting crap.
for station = 1:length(stationList)
    query = sprintf(['SELECT  msldepth, %s '                        ...
                     'FROM instrument.stationctd '                  ...
                     'WHERE time>=' '''%s''' ' and time<' '''%s'''  ...
                          ' and station=' '''%s'''                  ...
                     ' ORDER BY msldepth, time'],                   ...
                      varMap(variable),                             ... 
                      datestr(startDay, 'mmmm dd, yyyy'),           ...
                      datestr(endDay, 'mmmm dd, yyyy'),             ...
                      stationList{station});
    e = exec(dbConn,query);
    e = fetch(e);
    tempObsData.(stationList{station}) = e.data;
end

for station = 1:length(stationList)
    query = sprintf(['SELECT time '                                 ...
                     'FROM instrument.stationctd '                  ...
                     'WHERE time>=' '''%s''' ' and time<' '''%s'''  ...
                          ' and station=' '''%s'''                  ...
                     ' ORDER BY msldepth, time'],                   ...
                      datestr(startDay, 'mmmm dd, yyyy'),           ...
                      datestr(endDay, 'mmmm dd, yyyy'),             ...
                      stationList{station});
    e = exec(dbConn,query);
    e = fetch(e);
    obsTime.(stationList{station}) = e.data;
end

% Format observation data
%
% Change time to datenum, tempObsData to array, and combine in obsData which is
% a normal matrix, convert depth in centimeters to meters, and move different depths
% into structure  
% obsData.stationName.dep_n where n is the depth
%
for station = 1:length(stationList)
    times = NaN(size(obsTime.(stationList{station}),1),1);
    for row = 1:size(obsTime.(stationList{station}))
        times(row) = datenum(obsTime.(stationList{station}){row,1});
    end
    temp(1:size(times,1),1) = times(:,1);
    temp(1:size(times,1),2:3) = cell2mat(tempObsData.(stationList{station}));
    %temp(:,2) = temp(:,2)./100; Messes up struct names so names are in centimeters
    [b,m] = unique(temp(:,2));
    % If there are >1 depths seperate into different structs
    if size(b,1) > 1
        for depth = 1:size(b,1)
            if b(depth) == 99999
                continue;
            elseif depth == 1
                temp_2 = temp(1:m(depth),:);
            else
                temp_2 = temp(m(depth-1)+1:m(depth),:);
            end
            depStr = sprintf('dep_%i', b(depth));
            obsData.(stationList{station}).(depStr)(:,1) = temp_2(:,1);
            obsData.(stationList{station}).(depStr)(:,2) = temp_2(:,3);
        end
    % Dump data from one depth into struct
    else
        obsData.(stationList{station}).(sprintf('dep_%i',b(1)))(:,1) = temp(:,1);
        obsData.(stationList{station}).(sprintf('dep_%i',b(1)))(:,2) = temp(:,3);
    end
end


% Interpolate observation to model time step
% steptime in Julian date. 1 sec = 1/24/60/60, modelTs is in secs
interpStep = modelTs/24/60/60;
interpCut  = 3; 
startDay = addtodate(startDay, modelTs, 'second');
endDay

% Loop over each station and depth at each station
for station = 1:length(stationList)
    depths = fieldnames(obsData.(stationList{station}));
    for depth = 1:size(depths,1)
        % Remove NaNs to allow for interp to actually work
        idx = find(~isnan(obsData.(stationList{station}).(depths{depth})(:,2)));
        temp = obsData.(stationList{station}).(depths{depth})(idx,:); 
        temp = intdefine(temp(:,1), temp(:,2), interpStep, interpCut, startDay, endDay);
        obsInterp.(stationList{station}).(depths{depth}) = temp;
    end
end

