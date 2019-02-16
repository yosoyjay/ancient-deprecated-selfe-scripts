function [reg_kflow, reg_kwind, reg_vO2] = int_regions_NEM(kflow, kwind, vO2)
% [reg_kflow, reg_kwind, reg_vO2] = int_regions_NEM(int_kflow, int_kwind, int_vO2)
%
% This function integrates the temporally integrated values out of int_flow, etc. over 3 regions:
% 1. Columbia upstream of confluence with Willamette
% 2. Willamette upstream of confluence with Columbia
% 3. Columbia from Beaver Army upstream to confluence with Willamette
%
% Input:
%	kflow - kflow output from int_time_NEM.m [Integrated values over time] (cm/h) 
%	kwind - kwind " "
%	vO2 - vO2 " "
%
% Output:
%	reg_kflow - Hourly kflow values integrated over defined regions [unitless?] (cm/h * m^2) (days, hours, region)
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
reg_kflow = nan(nDays, nHours, nRegions);
reg_wind = nan(nDays, nHours, nRegions);
reg_vO2 = nan(nDays, nHours, nRegions);

% Calculate area for each element in each region
for region = 1:nRegions
	elemAreaInRegion(region,:) = regions(region).elems .* fg.ar'; % (region,nElements)
end

% Loop over days and hours and regions
for day = 1:nDays
%tic
	for hour = 1:nHours
	%tic
		% Convert node values to element values by averaging 3 nodes to element
		% Reshape data for given day and hour so it can be transposed
		% Data shaping tricks lifted from Mojy's habop
		shapedData = reshape(kflow(day,hour,:),1,size(kflow,3));
		dataElemNodes = shapedData(fg.e)';   % kflow -> (nNodes, nElements)
		kflowElems = nanmean(dataElemNodes); % kflowElems (1, nElements)

		shapedData = reshape(kwind(day,hour,:),1,size(kflow,3));
		dataElemNodes = shapedData(fg.e)';   % kwind -> (nNodes, nElements)
		kwindElems= nanmean(dataElemNodes);  % kwindElems (1, nElements)

		shapedData = reshape(vO2(day,hour,:),1,size(vO2,3));
		dataElemNodes = shapedData(fg.e)';   % kflow -> (nNodes, nElements)
		vO2Elems = nanmean(dataElemNodes);   % kflowElems (1, nElements)

		for region = 1:nRegions
		%tic
			% Integrate over each region
			% Elements in region * element area + filter out those not in regions (set to 0)
			% Multiply element area by kflowElem to get kflow for an element for an hour
			% kflowElem (cm/h) * elemAreaInRegion (m^2) HMMM...
			kflowElemArea = kflowElems .* elemAreaInRegion(region,:);	% (region,nea)       
			kwindElemArea = kwindElems .* elemAreaInRegion(region,:);
			vO2ElemArea = vO2Elems .* elemAreaInRegion(region,:);
			% Sum over entire region
			reg_kflow(day, hour, region) = nansum(kflowElemArea);
			reg_kwind(day, hour, region) = nansum(kwindElemArea);
			reg_vO2(day, hour, region) = nansum(vO2ElemArea); 
		%toc - ~0.005 seconds
		end
	%toc - ~ 0.05 seconds
	end
%toc - ~2 seconds per day
end


