function [kflow_min, kflow_max, kflow_mean] = avg_regions_NEM(kflow, kwind, vO2)
% [kflow_min, kwind_max, vO2_mean] = avg_regions_NEM(kflow, kwind, vO2)
%
% This function averages the temporally integrated values out of int_flow, etc. over 3 regions:
% 1. Columbia upstream of confluence with Willamette
% 2. Willamette upstream of confluence with Columbia
% 3. Columbia from Beaver Army upstream to confluence with Willamette
%
% Implmentation notes:
%	Nodal values for kflow, kwind, and vO2 are converted to elmental values to reuse code.
%   The difference between nodal and element values should be assessed.
%
% Input:
%	kflow - kflow output from int_time_NEM.m [Integrated values over time] (cm/h) 
%	kwind - kwind " "
%	vO2 - vO2 " "
%
% Output:
%	reg_kflow - Hourly kflow values integrated over defined regions (m/h * m^2) (days, hours, region)
% 	reg_kwind - " "
%	vO2 - vO2 - " "
%
% Hgrid and region paths are hard coded.
%
% lopezj - 3/18/2012
%

% Paths to add
addpath '/home/workspace/users/lopezj/scripts/matlab/regions/';

% Load grid information and river regions
grPath = '/home/workspace/users/lopezj/scripts/selfe_setup/db26/hgrid.gr3';
regPath = '/home/workspace/project/lopezj/projects/river/river_regions.dat';
[regions, fg] = read_regions(regPath, grPath);

% Consts used for dimensioning
nDays = size(kflow,1);
nHours = size(kflow,2);
nRegions = size(regions,2); 

% Allocate output arrays
kflow_min = zeros(nDays,nHours,nRegions);
kflow_max = zeros(nDays,nHours,nRegions);
kflow_mean = zeros(nDays,nHours,nRegions);

kwind_min = zeros(nDays,nHours,nRegions);
kwind_max = zeros(nDays,nHours,nRegions);
kwind_mean = zeros(nDays,nHours,nRegions);

vO2_min = zeros(nDays,nHours,nRegions);
vO2_max = zeros(nDays,nHours,nRegions);
vO2_mean = zeros(nDays,nHours,nRegions);

% Calculate area for each element in each region to filter elements later.
% If element is not in a region it has an area of 0.
% If element is in region give it a value of 1.
%for region = 1:nRegions
%	elemAreaInRegion(region,:) = regions(region).elems .* fg.ar'; % (region,nElements)
%end
%elemAreaInRegion(elemAreaInRegion~=0) = 1;

% Loop over days and hours and regions
for day = 1:nDays
%tic
	for hour = 1:nHours
	%tic
		% Convert node values to element values by averaging 3 nodes to element
		% Reshape data for given day and hour so it can be transposed
		% Data shaping tricks lifted from Mojy's habop
		%shapedData = reshape(kflow(day,hour,:),1,size(kflow,3));
		%dataElemNodes = shapedData(fg.e)';   % kflow -> (nNodes, nElements)
		%kflowElems = nanmean(dataElemNodes); % kflowElems (1, nElements)

		%shapedData = reshape(kwind(day,hour,:),1,size(kflow,3));
		%dataElemNodes = shapedData(fg.e)';   % kwind -> (nNodes, nElements)
		%kwindElems= nanmean(dataElemNodes);  % kwindElems (1, nElements)

		%shapedData = reshape(vO2(day,hour,:),1,size(vO2,3));
		%dataElemNodes = shapedData(fg.e)';   % kflow -> (nNodes, nElements)
		%vO2Elems = nanmean(dataElemNodes);   % kflowElems (1, nElements)

		for region = 1:nRegions
		%tic
			% Average over each region
			% 1. Id each node in the region
			% 2. Filter shaped data so only nodes inside region are not zero
			% 3. Change all zeros to NaNs.
			inRegion = inpolygon(fg.x                       , fg.y, ...
							     regions(region).points(:,1), regions(region).points(:,2));

			shapedData = reshape(kflow(day,hour,:),1,size(kflow,3));   % (1,np)
		    shapedData = shapedData(inRegion);
		    kflow_min(region) = min(shapedData);
			kflow_max(region) = max(shapedData);
			kflow_mean(region) = nanmean(shapedData);

			shapedData = reshape(kwind(day,hour,:),1,size(kflow,3));   % (1,np)
		    shapedData = shapedData(inRegion);
		    kflow_min(region) = min(shapedData);
			kflow_max(region) = max(shapedData);
			kflow_mean(region) = nanmean(shapedData);

			shapedData = reshape(vO2(day,hour,:),1,size(kflow,3));   % (1,np)
		    shapedData = shapedData(inRegion);
		    kflow_min(region) = min(shapedData);
			kflow_max(region) = max(shapedData);
			kflow_mean(region) = nanmean(shapedData);
		%toc - ~0.005 seconds
		end
	%toc - ~ 0.05 seconds
	end
%toc - ~2 seconds per day
end


