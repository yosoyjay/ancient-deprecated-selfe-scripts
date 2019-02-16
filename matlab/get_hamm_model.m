function [obs, model] = get_hamm_model(startDate, endDate, modelPath)
% Extracts observation and model results and returns it with time aligned.
%
% Input:
%   startDate - Datenum
%   endDate   - Datenum
%   modelPath - Path to run directory with model results
% Output:
%   obs       - Observation data from Hammond tide gauge
%   model     - Model results for modelPath run
%
% lopezj - 11/21/2011

% Prep work of dates to "day of year"
sd = datevec(startDate);
ed = datevec(endDate);
start = date2jd(sd(1), sd(2), sd(3), sd(4), sd(5), sd(6));
stop   = date2jd(ed(1), ed(2), ed(3), ed(4), ed(5), ed(6)); 

% Load Hammond data
load '/home/workspace/project/lopezj/data/hammond/hammond_Aug_2010';
time = data.values(:,1);
elev = data.values(:,2);
hamm = [-123.95, 46.2];

% Set up dates for interplotion and model data extraction 
date(1) = startDate;
i = 1;
while date(i) <= endDate
	i = i+1;
	date(i) = addtodate(date(i-1), 15, 'minute');
end 

% Find dates of interest then interpolate observation values in time to 15 min 
interpStep = 15/24/60;
interpCut  = 3;					 
times = find(time > startDate & time < endDate);
interpElev = intdefine(time(times(1):times(end)), elev(times(1):times(end)), ...
			           interpStep, interpCut, startDate, endDate);
interpTime = interpElev(:,1);

% Only get model data that I want lined up with depths lined up from above.  This takes forever
% depending on the size of the data.
for i=1:size(interpTime,1)
	modData = modext(interpTime(i,1), hamm(1,2), hamm(1,1), interpElev(i,2), 'Elev', modelPath);
end

model = modData;
obs   = interpElev;

