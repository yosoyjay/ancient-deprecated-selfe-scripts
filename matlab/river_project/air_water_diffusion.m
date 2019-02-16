function [k_flow, k_wind, vO2, eta] = air_water_diffusion(runPath, startDay, endDay)
% [k_flow, k_wind, vO2] = air_water_diffusion(runPath, startDay, endDay)
%
% This function calculates the air water diffusion based on model results
% of water velocity, river depth, temperature, and wind speed.
%
% Input:
%	runPath - Directory to model run results
%	startDay - Day to start calculating. This is integer of model run day.
%	endDay - Last day to caclculate.  This is an interger of model run day.
%
% Output:
%	k_flow(day, ts, np) - Water component of the air-water diffusion (cm/h)
%	k_wind(day, ts, np - Wind component of the air-water diffusion (cm/h)
%	vO2(day, ts, np) - The sum of k_flow and a_flow. (cm/h)
%
% lopezj - 03/11/2012
%
tic
% M-elio library used to extract model results.
addpath /usr/local/cmop/matlab/cmop/m-elio/;
% calc_mole_diff is here
addpath /home/workspace/users/lopezj/scripts/matlab/river_project/;
addpath /home/workspace/users/lopezj/scripts/matlab/;
% Extract number of time steps per file
addpath /home/workspace/users/lopezj/scripts/matlab/model_extraction/

% Assume there are 96 time steps in output
% timeSteps = [1:2:96]; 

% Deal with h- and v-grid
hgridPath = sprintf('%s/hgrid.gr3', runPath);
gr.hgrid = gr_readHGrid(hgridPath);


% Loop over all days and calculate k_flow, a_flow, and vO2 at every node
for day = startDay:endDay
	% Load headers for model output files
	hElev = sz_readHeader([runPath '/outputs/' num2str(day) '_elev.61']);
	hDahv = sz_readHeader([runPath '/outputs/' num2str(day) '_dahv.62']);
	hTemp = sz_readHeader([runPath '/outputs/' num2str(day) '_temp.63']);
	hWind = sz_readHeader([runPath '/outputs/' num2str(day) '_wind.62']);
	gr.vgrid = hElev.vgrid;

	% Get time steps and number of steps information
	dt = hElev.dt;
	nSteps = hElev.nSteps;

	% Loop over time steps for each day every half hour for now.
	for timeStep = 1:2:nSteps
		% Read the time steps for the current day
		[dElevs ts] = sz_readTimeStep(hElev,timeStep);		
		[dTemps ts] = sz_readTimeStep(hTemp,timeStep);		
		[dHvels ts] = sz_readTimeStep(hDahv,timeStep);		
		[dWinds ts] = sz_readTimeStep(hWind,timeStep);		

		% Map sz levels to depths if required (3d vars only)
		elevs = dElevs; 
		temps = map_sz2hts(hTemp, dTemps(:));
		hvels = dHvels; 
		winds = dWinds; 

		% Construct vertical grid - Not necessary because 3d vars are averaged
		% sz = sz_computeZlevels(gr.hgrid.depth, elevs, gr.vgrid);	

		% Calculate depth at every node... if < -9000 it's dry.
		depths = gr.hgrid.depth + elevs;
		depths(depths<-9000) = NaN;

		% Find average temperature for the water column	and calculate molecular diffusion
		temps(temps==-99)=NaN;
		t_mean = nanmean(temps')';
		[diff,eta] = calc_mole_diff(t_mean);

		% Magnitude of velocity vectors
		velocity = sqrt(hvels(:,1,:).*hvels(:,1,:)+hvels(:,2,:).*hvels(:,2,:));
		winds_mag = sqrt(winds(:,1,:).*winds(:,1,:)+winds(:,2,:).*winds(:,2,:));

		% Unit conversions
		% Depths in centimeters
		depths = depths.*100;
		% Velocity in cm s^-1
		velocity = velocity.*100;
		
		% Calculate k_flow
		k_flow(day, timeStep, :) = sqrt(velocity.*diff./depths);

		% Calculate air_flow using Schmidt number
		sc_0 = 1800.6 - 120.1.*t_mean + 3.7818.*t_mean.*t_mean - 0.047608.*t_mean.*t_mean.*t_mean;
		salt = 0; % Salt is assumed to be zero because we are only looking upstream of Beaver Army
				  % Can easily be manipulated to load salt values and calculated for estuary
		sc = sc_0.*(1+3.14*10^-3*salt);
		k_wind(day, timeStep, :) = 0.31.*winds_mag.*winds_mag.*(sc./660).^(-0.5);	

		% Calculate vO2
		vO2(day, timeStep, :) = k_flow(day, timeStep, :) + k_wind(day, timeStep, :);
	end
end
toc
