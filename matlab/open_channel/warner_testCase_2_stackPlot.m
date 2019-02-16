function [] = windMixingTransect(runPath, movie)
% This function creates a plot for a transects similar to those found in Warner, but using
% filled contours.
%
% Input:     
%   runPath - Path to run directory, assumes outputs are already combined in outputs dir
%           - Assumes *.bp file named burchard.bp in directory
%			- Time step # (Not time) to plot data? 
%   movie   - 1 - Makes a avi, 0 - No avi
%
% Output:
%
% lopezj - 12/08/11
%
tic;

% Constants etc.
%addpath('/usr/local/cmop/matlab/cmop/m-elio');
addpath('/home/workspace/users/lopezj/svn/modeling/matlab/m-elio/trunk');
addpath('/home/workspace/users/lopezj/scripts/matlab');
BP_FILE = 'oc_length.bp';
MOVIE_FILE = 'wind_mixing.avi';

% Make a plots folder if it doesn't exists for the output 
plots = sprintf('%s/plots', runPath);
if (~exist(plots,'file'))
	mkdir(plots);
end

% Build points file
bpPath = sprintf('%s/%s', runPath, BP_FILE);
tr.hgrid = gr_readHGrid(bpPath);

% Create animation 
if movie == 1
	moviePath = sprintf('%s/%s', plots, MOVIE_FILE);
	movieH = avifile(moviePath);
end

% Get information about the number of days in run from param.in
% Assumes it is on line 151
paramPath = sprintf('%s/param.in', runPath);
param = fopen(paramPath, 'r');
C = textscan(param,'%s','delimiter','\n');
dayRawStr = C{1,1}(151);
dayRawStr = regexp(dayRawStr,' ','split');
days = str2num(cell2mat(dayRawStr{1,1}(3)));

% Axis limits for all variables
axisDimsTdff = [0 0.025];
axisDimsVdff = [0 0.025];
axisDimsHvel = [-1 1];
axisDimsKine = [0 0.0005];
axisDimsSed  = [0 1];
axisDimsSalt = [0 30];

% X Axis ticks and labels based on build points file
% xRange set to show 1/n of the domain!
xMin = roundsd(min(tr.hgrid.x),3);
xMax = roundsd(max(tr.hgrid.x),3);
xRange = floor((xMax-xMin)/4);
xRangeTicks = floor(xRange/1000);
xTicks = xRangeTicks/4;

% Axis for transect x & z and label
axisDims = [xMin xMax -30 3];
aTTicks   = xMin:xRange:xMax;
aTLabel   = 0:xRangeTicks:xRangeTicks*4;

for day=1:10
	% Deal with h- and v-grid
	% Note that sz layers must be extracted from binary output file headers
	% Open up files for all variables I'm looking at in open channel
	hgridPath = sprintf('%s/hgrid.gr3', runPath);
	gr.hgrid = gr_readHGrid(hgridPath);
	hvelPath = sprintf('%s/outputs/%d_hvel.64', runPath, day);
	hHvel = sz_readHeader(hvelPath);
	kinePath = sprintf('%s/outputs/%d_kine.63', runPath, day);
	hKine = sz_readHeader(kinePath);
	tdffPath = sprintf('%s/outputs/%d_tdff.63', runPath, day);
	hTdff = sz_readHeader(tdffPath);
	vdffPath = sprintf('%s/outputs/%d_vdff.63', runPath, day);
	hVdff = sz_readHeader(vdffPath);
	elevPath = sprintf('%s/outputs/%d_elev.61', runPath, day);
	hElev = sz_readHeader(elevPath);
	saltPath = sprintf('%s/outputs/%d_salt.63', runPath, day);
	hSalt = sz_readHeader(saltPath);
	% sedPath = sprintf('%s/outputs/%d_trcr_1.63', runPath, day);
	% hSed = sz_readHeader(sedPath);
	gr.vgrid = hKine.vgrid;

	% Compute transect
	[ob]= ob_ini_fromTrasect(gr, bpPath);
	trLen = cumsum(sqrt((ob.xy.x(2:end)-ob.xy.x(1:end-1)).^2 + ...
						(ob.xy.y(2:end)-ob.xy.y(1:end-1)).^2));
	trLen = [0; trLen];

	% Create a plot for every time step
	% Get time step and number of step information
	dt = hHvel.dt;
	nSteps = hHvel.nSteps;
  dt
  nSteps
	for i=1:nSteps
		% Read timestep of data
		[hvelData varTs] = sz_readTimeStep(hHvel,i);
		[kineData varTs] = sz_readTimeStep(hKine,i);
		[tdffData varTs] = sz_readTimeStep(hTdff,i);
		[vdffData varTs] = sz_readTimeStep(hVdff,i);
		[saltData varTs] = sz_readTimeStep(hSalt,i);
		% [sedData varTs]  = sz_readTimeStep(hSed,i);

		% Map sz levels to depths
		dHvel = map_sz2hts(hHvel, hvelData(:,1));
		dKine = map_sz2hts(hKine, kineData);
		dTdff = map_sz2hts(hTdff, tdffData);
		dVdff = map_sz2hts(hVdff, vdffData);
		dSalt = map_sz2hts(hSalt, saltData);
		% dSed  = map_sz2hts(hSed, sedData);

		% Do this and you just end up with the data that you want in transect
		% ob is the object that has the transect information
		hvel = ob.xy.H*double(dHvel);
		kine = ob.xy.H*double(dKine);
		tdff = ob.xy.H*double(dTdff);
		vdff = ob.xy.H*double(dVdff);
		salt = ob.xy.H*double(dSalt);
		% sed  = ob.xy.H*double(dSed);

		% Construct vertical grid
		% Read timestep of elevations
		elevData = sz_readTimeStep(hElev,i);
		elev     = ob.xy.H*double(elevData);
		depths   = ob.xy.H*gr.hgrid.depth;
		sz=sz_computeZlevels(depths,elev,gr.vgrid);

		% Copy x values so every node vertical level has an x value
		x = repmat(ob.xy.x, 1, size(sz,2)); 

		% Plot data and analytical solution
		%screen_sz = get(0,'ScreenSize');
		%figH = figure('Position',[ 1 screen_sz(4)/2 1200 800]);
		figH = figure('visible', 'off','color','w');

		% Velocity
		h(1) = subplot(5,1,1);
		sp_pos = get(h(1), 'position'); hold on;
		contourf(x,sz,-hvel,'linecolor','none');
%		colormap(bluewhitered);
		colorbar;
		caxis(axisDimsHvel);
		title('Velocity m/s');
		%legH = legend('Model', 'Analytical', 'Location', 'NorthEast');
		ylabel('Depth m');
    plot(x(:,1:3), sz(:,1:3), 'color', [0.7 0.7 0.7], 'linewidth', 0.1);
    axis(axisDims);
    set(gca,'xtick',aTTicks);
    set(gca,'xticklabel',aTLabel);


		% Salt
		h(2) = subplot(5,1,2);
		sp_pos = get(h(2), 'position'); hold on;
		contourf(x,sz,salt,'linecolor','none');
		colorbar;
		caxis([axisDimsSalt]);
    plot(x(:,1:3), sz(:,1:3), 'color', [0.7 0.7 0.7], 'linewidth', 0.1);
		title('Salt psu');
    axis(axisDims);
    set(gca,'xtick',aTTicks);
    set(gca,'xticklabel',aTLabel);

		% TKE
		h(3) = subplot(5,1,3);
		sp_pos = get(h(3), 'position'); hold on;
		contourf(x,sz,kine,'linecolor','none');
		axis(axisDims);
    plot(x(:,1:3), sz(:,1:3), 'color', [0.7 0.7 0.7], 'linewidth', 0.1);
		title('TKE m^2/s^2');
		ylabel('Depth m');
		%colormap(jet(1024));
		%colorbar_log([0.0001 0.1]);
		colorbar;
		caxis(axisDimsKine);
    axis(axisDims);
    set(gca,'xtick',aTTicks);
    set(gca,'xticklabel',aTLabel);

		% Plot levels for reference
		%plot(x(:,end), sz(:,end), 'b', 'linewidth', 0.1, 'linesmoothing', 'on'); hold on;
		%plot(x, sz, 'color', 'w','linewidth', 0.1, 'linesmoothing', 'on'); hold on;
		%axis(axisDims);

		% Eddy diffusivity 
		h(4) = subplot(5,1,4);
		sp_pos = get(h(4), 'position'); hold on;
		contourf(x,sz,tdff,'linecolor','none');
		title('Diffusivity m^2/s');
		axis(axisDims);
    plot(x(:,1:3), sz(:,1:3), 'color', [0.7 0.7 0.7], 'linewidth', 0.1);
		%colormap(jet(1024));
		%colorbar_log([0.01 0.1]);
		colorbar;
    axis(axisDims);
    set(gca,'xtick',aTTicks);
    set(gca,'xticklabel',aTLabel);

		% Eddy viscosity 
		h(5) = subplot(5,1,5);
		sp_pos = get(h(5), 'position'); hold on;
		contourf(x,sz,vdff,'linecolor','none');
		title('Viscosity m^2/s');
		axis(axisDims);
    plot(x(:,1:3), sz(:,1:3), 'color', [0.7 0.7 0.7], 'linewidth', 0.1);
		%colormap(jet(1024));
		%colorbar_log([0.01 0.1]);
		colorbar;
    axis(axisDims);
    set(gca,'xtick',aTTicks);
    set(gca,'xticklabel',aTLabel);

		% Finally add the time and text for colorbar
		% tah - text axis handle
		% ts  - timeString
		% th  - text handle
		secs = (day-1)*dt*nSteps + i*dt;
		tideCycle = secs/86400;
		tah = axes('position', [0 0 1 1]);
		ts  = sprintf('%2.2f %s', tideCycle, 'days');
		th  = text(0.5, 0.98, ts, 'FontSize', 12);
		set(th, 'HorizontalAlignment', 'center');
		th  = text(0.5, 0.06, 'Distance from ocean boundary km', 'FontSize', 12);
		set(th, 'HorizontalAlignment', 'center');
		set(gca,'Visible', 'off');

		% Get ready for saves
		fn = sprintf('%s/%s-%06d.png', plots, 'wind_mixed', secs); 
		iw = 2048;
		ih = 1600;
		set(gcf,'PaperUnits','inches','PaperPosition',[0 0 iw/100.0 ih/100.0])
		print('-dpng', fn, '-r100');
		if movie == 1
			movieH = addframe(movieH, figH);
		end
		close(figH);
	end % Loop over i = 1:timeSteps - Time steps per day
end	% Loop for day = 1:rnday - The output days
if movie == 1
    movieH = close(movieH);
end
