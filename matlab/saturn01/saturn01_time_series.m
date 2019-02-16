function [] = station_time_series(runPath, bpFile, startDay, endDay, variable)
% Creates a time series plot of a single point for the given variable
%
% Input:
%	runPath - Path to the run directory, assumes outputs are already combined in outputs dir
%	bpFile  - build point file of the station
%	startDay - Day number of run (File number 1_salt, 2_elev, etc...)
%	endDay   - Day number of run (File number 1_salt, 2_elev, etc...)
%	variable - Variable to plot ('salt.63', 'trcr_1.63', 'temp.63', etc..)
%
% lopezj
%

% Constants etc.
addpath('/usr/local/cmop/matlab/cmop/m-elio');
addpath('/home/workspace/users/lopezj/scripts/matlab');

% Make a plots folder if it doesn't exists for the output 
plots = sprintf('%s/plots', runPath);
if (~exist(plots,'file'))
	mkdir(plots);
end

% Axis limits for all variables
axisDimsTdff = [10e-6 0];
axisDimsVdff = [10e-6 0];
axisDimsHvel = [-2 2];
axisDimsKine = [0 0.15];
axisDimsSed  = [0 0.1];
axisDimsSalt = [0 35];

% Get information about the run
runInfo = get_run_info(runPath);

% Load up build points file
tr.hgrid = gr_readHGrid(bpFile);

% Extract model data from all the days
for day = startDay:endDay
	% Deal with h- and v-grid
	% Note that sz layers must be extracted from binary output file headers
	% Open up files for the variable
	hgridPath = sprintf('%s/hgrid.gr3', runPath);
	gr.hgrid = gr_readHGrid(hgridPath);
	elevPath = sprintf('%s/outputs/%d_elev.61', runPath, day);
	hElev = sz_readHeader(elevPath);
	varPath = sprintf('%s/outputs/%d_%s', runPath, day, variable);
	hVar = sz_readHeader(varPath);
	gr.vgrid = hVar.vgrid;

	% Compute transect
	[ob]= ob_ini_fromTrasect(gr, bpFile);
	trLen = cumsum(sqrt((ob.xy.x(2:end)-ob.xy.x(1:end-1)).^2 + ...
						(ob.xy.y(2:end)-ob.xy.y(1:end-1)).^2));
	trLen = [0; trLen];

	% Collect data for every time
	[varData varTs] = sz_readTimeStep(hVar, 1:runInfo.nSteps); 
	
	% Map sz levels to depths
	dVar = map_sz2hts(hVar, varData);

	% Limits data to just that in the transect
	var = ob.xy.H*double(dVar);

	% Construct vertical grid
	elevData = sz_readTimeStep(hElev, 1:runInfo.nSteps);
	elev     = ob.xy.H*double(elevData);
	depths   = ob.xy.H*gr.hgrid.depth;
	sz = sz_computeZLevels(depths, elev, gr.vgrid);

	% Copy x values so every nodes vertical level has an x value
	x = repmat(ob.xy.x, 1, size(sz,2));
end

% Create the plot!
figH = figure('visible', 'off', 'color', 'w');
contourf(x, sz, var, 'linecolor', 'none', 'levelstep', 2);
box on;
colorbar;
caxis([0 32]);
title('Salt psu');
h = get(gca,'title');
set(gca,'fontsize',14);
fn = sprintf('%s/%s.png', plots, 'n_channel_test'); 
iw = 1024;
ih = 800;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 iw/100.0 ih/100.0])
print('-dpng', fn, '-r100');


