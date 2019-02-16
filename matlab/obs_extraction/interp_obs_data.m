function [obsInterp] = interp_obs_data(obsData, startTime, endTime, modelTs)
% [obsInterp] = interp_obs_data(obsData, startTime, endTime, modelTs)
% Interploates observation data to model time step for skill assessment
%
% Input:
%	odsData - Column 1 time in datenum, column 2 data 
%	startTime - Time when model data starts in datenum format
%	endTime - Time when model data ends in datenum
%	modelTS - Output TS
%	
% Output:
%	obsInterp - Observation data interploated to modelTimeStep
%
%
interpStep = modelTs/24/60/60;
interpCut  = 3; 

% Remove NaNs to allow for interp to actually work
idx = find(~isnan(obsData(:,2)));
temp = obsData(idx,:); 

% Remove potential for duplicate times that causes intdefine to die
[b,m] = unique(temp(:,1));
temp = temp(m,:);
temp = intdefine(temp(:,1), temp(:,2), interpStep, interpCut, startTime, endTime);

%
obsInterp = temp;


