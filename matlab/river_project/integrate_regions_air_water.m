function [reg_kflow, reg_kwind, reg_vO2] = integrae_regions_air_water(int_kflow, int_kwind, int_vO2)
% [reg_kflow, reg_kwind, reg_vO2] = integrae_regions_air_water(int_kflow, int_kwind, int_vO2)
%
% This function integrates the temporally integrated values out of int_flow, etc. over 3 regions:
% 1. Columbia upstream of confluence with Willamette
% 2. Willamette upstream of confluence with Columbia
% 3. Columbia from Beaver Army upstream to confluence with Willamette
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
nDays = size(int_flow,1);
nHours = size(int_flow,2);
nRegions = size(regions,2); 

% Allocate output arrays
reg_kflow = nan(nDays, nHours, nRegions);
reg_wind = nan(nDays, nHours, nRegions);
reg_vO2 = nan(nDays, nHours, nRegions);

% Loop over days and hours and regions
for day = 1:nDays
	for hour = 1:nHours
		for regions = 1:nRegions
			% Convert node values to element values by averaging 3 nodes to element
			kflow_at_nodes_in_elems = kflow(day,hour,fg.e);
			kflow_elems = nanmean(kflow_at_nodes_in_elems');
			% Integrate over each region
			% Calculate k_flow at element (average of three nodes)
			% sum(kflow*100*area of element)
			idx = regions(region).elems .* kflow_elems; 
			reg_kflow(day, hour, :) = 
			reg_wind(day, hour, :) = 
			reg_vO2(day, hour, :) = 
		end
	end
end


