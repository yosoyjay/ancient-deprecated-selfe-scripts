function [sz, sed] = settling_plot(runPath, startDay, endDay, plots) 
% [sz, sed] = settling_plot(runPath, startDay, endDay) 
%
% This function extracts sediment data, creates a plot, and returns
% the vertical levels(sz) and [sediment]
%
% Input:
%	runPath - Directory to model run results
%	startDay - Day to start calculating. This is integer of model run day.
%	endDay - Last day to caclculate.  This is an interger of model run day.
%	plots - 1==yes, 0==no
%
% Output:
%   sz - Vertical level information	
%	sed - Sediment concentrations
%
% lopezj - 04/05/2012
%

% M-elio library used to extract model results.
addpath /usr/local/cmop/matlab/cmop/m-elio/;

% Extract number of time steps per file
addpath /home/workspace/users/lopezj/scripts/matlab/model_extraction/;
addpath /home/workspace/users/lopezj/scripts/matlab/;

% Read every skip output time steps
skip = 1;

% Deal with h-grid 
hgridPath = sprintf('%s/hgrid.gr3', runPath);
gr.hgrid = gr_readHGrid(hgridPath);

% Set up plots
if plots == 1
	mkdir('plots');
end


% Loop over all days and extract data 
for day = startDay:endDay
	% Load headers for model output files
	hElev = sz_readHeader([runPath '/outputs/' num2str(day) '_elev.61']);
	hSed1 = sz_readHeader([runPath '/outputs/' num2str(day) '_trcr_1.63']);
	gr.vgrid = hElev.vgrid;

	% Get time steps and number of steps information
	dt = hElev.dt;
	nSteps = hElev.nSteps;

	% Define output times steps to read
	timeStep = 1:skip:nSteps;

	% Read the time steps for the current day
	% Returned data shape (nodes, depths, timeStep)
	[dElevs ts] = sz_readTimeStep(hElev,timeStep);		
	[dSeds1 ts] = sz_readTimeStep(hSed1,timeStep);		

	% Map sz levels to depths if required (3d vars only)
	elevs = dElevs; 
	sed = map_sz2hts_mat(hSed1, dSeds1);

	% Construct vertical grid - Not necessary because 3d vars are averaged
	sz = sz_computeZlevels(gr.hgrid.depth, reshape(elevs(:,1,1),size(elevs,1),1), gr.vgrid);	

	% Calculate depth at every node... if = -9999 it's dry.
	% Repeat depths to match shape of elevs mat.
	depths = repmat(gr.hgrid.depth,[1 size(elevs,2) size(elevs,3)]) + elevs;
	depths(depths==-9999) = NaN;

	% Create plots
	for ts = 1:size(elevs,3)
		if ts == 1
			color = [0 0 0];
		else
			color = [0.6 0.6 0.6];
		end
		plot(sed(1,:,ts),sz(1,:),'color',color);	hold on;
		axis([0, 0.125,-1,0]);
		box on;
	end
end
%toc
