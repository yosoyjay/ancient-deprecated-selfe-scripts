function [obsData] = get_obs_data(query);
% [obsData] = get_obs_data(query);
% Convenience function to do very simple queries 
%
% Input:
%	vars  - Variables for SELECT (cell array)
%	table - Table for FROM (string)
%	startTime - Time to start query (datenum)
%	endTime - Time to end query (datenum)
%	stat  - Station to get data from (string, if applicable)
% 	order - ORDER BY 
% Output:
%	obsData(n,1+size(vars)) - Array of data (time, data)
%
% lopezj - 01/25/2012
%

% Connect to database
dbConn = connect_to_db();

% Seperate time and data in two different queries to minimize messing with Matlab
% formatting crap.
for station = 1:length(stationList)
    query = sprintf(['SELECT  %s '                                  ...
                     'FROM %s '                                     ...
                     'WHERE timezone(''PST'',time)>=' '''%s'''      ...
					      ' and timezone(''PST'',time)<' '''%s'''   ...
                          ' and station=' '''%s'''                  ...
                     ' ORDER BY %s'],                               ...
					  varMap(variable).var,							...
                      varMap(variable).table,                       ... 
                      datestr(startDay, 'mmmm dd, yyyy'),           ...
                      datestr(endDay, 'mmmm dd, yyyy'),             ...
                      stationList{station},                         ...
					  varMap(variable).order);
    e = exec(dbConn,query);
    e = fetch(e);
    tempObsData.(stationList{station}) = e.data;
end

for station = 1:length(stationList)
    query = sprintf(['SELECT timezone(''PST'',time) '                ...
                     'FROM %s '                                     ...
                     'WHERE timezone(''PST'',time)>=' '''%s'''      ... 
					      ' and timezone(''PST'',time)<' '''%s'''   ...
                          ' and station=' '''%s'''                  ...
                     ' ORDER BY %s'],                               ...
					  varMap(variable).table,                       ...
                      datestr(startDay, 'mmmm dd, yyyy'),           ...
                      datestr(endDay, 'mmmm dd, yyyy'),             ...
                      stationList{station},                         ...
					  varMap(variable).order);
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
	if variable == 'elev.61'
    	depStr = sprintf('dep_%i', 0);
        obsData.(stationList{station}).(depStr)(:,1) = temp(:,1);
        obsData.(stationList{station}).(depStr)(:,2) = temp(:,3);
	else	
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
end
 Connect to database
dbConn = connect_to_db();

% Seperate time and data in two different queries to minimize messing with Matlab
% formatting crap.
for station = 1:length(stationList)
    query = sprintf(['SELECT  %s '                                  ...
                     'FROM %s '                                     ...
                     'WHERE timezone(''PST'',time)>=' '''%s'''      ...
					      ' and timezone(''PST'',time)<' '''%s'''   ...
                          ' and station=' '''%s'''                  ...
                     ' ORDER BY %s'],                               ...
					  varMap(variable).var,							...
                      varMap(variable).table,                       ... 
                      datestr(startDay, 'mmmm dd, yyyy'),           ...
                      datestr(endDay, 'mmmm dd, yyyy'),             ...
                      stationList{station},                         ...
					  varMap(variable).order);
    e = exec(dbConn,query);
    e = fetch(e);
    tempObsData.(stationList{station}) = e.data;
end

for station = 1:length(stationList)
    query = sprintf(['SELECT timezone(''PST'',time) '                ...
                     'FROM %s '                                     ...
                     'WHERE timezone(''PST'',time)>=' '''%s'''      ... 
					      ' and timezone(''PST'',time)<' '''%s'''   ...
                          ' and station=' '''%s'''                  ...
                     ' ORDER BY %s'],                               ...
					  varMap(variable).table,                       ...
                      datestr(startDay, 'mmmm dd, yyyy'),           ...
                      datestr(endDay, 'mmmm dd, yyyy'),             ...
                      stationList{station},                         ...
					  varMap(variable).order);
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
	if variable == 'elev.61'
    	depStr = sprintf('dep_%i', 0);
        obsData.(stationList{station}).(depStr)(:,1) = temp(:,1);
        obsData.(stationList{station}).(depStr)(:,2) = temp(:,3);
	else	
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
end

