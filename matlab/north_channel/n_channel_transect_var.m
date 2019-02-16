function [] = n_channel_transect(runPath, bpFile, startDay, endDay, var)
% [] = n_channel_transect(runPath, bpFile, startDay, endDay, var) 
% This function creates a plot for a transect of the North Channel or anywhere as long as
% you give it a *.bp file.
%
% Made to look at nSedClasses and salinity in the north channel.
%
% Input:     
%   runPath  - Path to run directory, assumes outputs are already combined in outputs dir
%            - Assumes *.bp file named n_channel.bp in directory
%   bpFile   - *.bp file relative to runPath
%   startDay - Day number of run (File number 1_hvel, 2_hvel, etc...)
%   endDay   - Day number of run (File number 6_hvel, 12_hvel, etc...)
%		var 	 - Output variable in form of kine.63, salt.63, temp.63 etc.
%
% Output:
%
% lopezj - 01/05/12
% lopezj - 04/10/12 - Changed to handle variable number of sediment classes
% lopezj - 06/12/12 - Changed to handle any variable and change file naming scheme
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
axisDimsKine = [0 0.01];
axisDimsSed  = [0 0.1];
axisDimsSalt = [0 32];
axisDims = axisDimsKine;

% Build points file
bpPath = sprintf('%s/%s', runPath, bpFile);
tr.hgrid = gr_readHGrid(bpPath);

% X Axis ticks and labels based on build points file
% xRange set to show 1/n of the domain!
xMin = roundsd(min(tr.hgrid.x),3);
xMax = roundsd(max(tr.hgrid.x),3);
xRange = floor((xMax-xMin)/4);
xRangeTicks = floor(xRange/1000);
xTicks = xRangeTicks/4;

% Axis for transect x & z and label
axisTrans = [xMin xMax -30 3];
%axisTrans = [xMin xMax -15 3];
aTTicks   = xMin:xRange:xMax;
aTLabel   = 0:xRangeTicks:xRangeTicks*4;

% Get information about the number of days in run from param.in
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% IMPORTANT Assumes it is on line 151		IMPORTANT
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
assumedLine = 151;		% If the line number of rnday changes script will day --- change here
paramPath = sprintf('%s/param.in', runPath);
param = fopen(paramPath, 'r');
C = textscan(param,'%s','delimiter','\n');
dayRawStr = C{1,1}(assumedLine);
dayRawStr = regexp(dayRawStr,' ','split');
days = str2num(cell2mat(dayRawStr{1,1}(3)));
fclose(param);

% Get information about the start date from bctides - lifted from Grant's stuff
%bctidesPath = sprintf('%s/bctides.in', '../');
bctidesPath = sprintf('%s/bctides.in', runPath);
fid = fopen(bctidesPath, 'r');
dstr = fgets(fid);
startDate = datenum(dstr(1:end-4),'mm/dd/yyyy HH:MM:SS');
fclose(fid);
                                                                                                        
for day=startDay:endDay
	tic
	% Deal with h- and v-grid
	% Note that sz layers must be extracted from binary output file headers
	% Open up files for all variables I'm looking at in open channel
	hgridPath = sprintf('%s/hgrid.gr3', runPath);
	gr.hgrid = gr_readHGrid(hgridPath);
	elevPath = sprintf('%s/outputs/%d_elev.61', runPath, day);
	hElev = sz_readHeader(elevPath);
	varPath = sprintf('%s/outputs/%d_%s', runPath, day,var);
	hVar = sz_readHeader(varPath);
	gr.vgrid = hVar.vgrid;

	% Compute transect
	[ob]= ob_ini_fromTrasect(gr, bpPath);
	trLen = cumsum(sqrt((ob.xy.x(2:end)-ob.xy.x(1:end-1)).^2 + ...
						(ob.xy.y(2:end)-ob.xy.y(1:end-1)).^2));
	trLen = [0; trLen];

	% Get time step and number of step information
	dt = hVar.dt;
	nSteps = hVar.nSteps;

%   Get all data for the day
%   [saltData varTs] = sz_readTimeStep(hVar,1:nSteps);
%   [elevData varTs] = sz_readTimeStep(hElev,1:nSteps);
%
%   Map sz levels to depths for 3d vars
%   dVar = map_sz2hts_mat(hVar, saltData)

	% Create a plot for every n time steps
	for i=1:nSteps
%       % Get data for day 
		%[hvelData varTs] = sz_readTimeStep(hHvel,i);
		[varData varTs] = sz_readTimeStep(hVar,i);
		%[sedData varTs]  = sz_readTimeStep(hSed,i);

		% Map sz levels to depths
		dVar = map_sz2hts(hVar, varData);
		%dSed  = map_sz2hts(hSed, sedData);
		%dHvel = map_sz2hts(hHvel, hvelData(:,1));

		% Do this and you just end up with the data that you want in transect
		%hvel = ob.xy.H*double(dHvel);
		var = ob.xy.H*double(dVar);
		%sed  = ob.xy.H*double(dSed);

		% Construct vertical grid
		% Read timestep of elevations
		elevData = sz_readTimeStep(hElev,i);
		elev     = ob.xy.H*double(elevData);
		depths   = ob.xy.H*gr.hgrid.depth;
		sz=sz_computeZlevels(depths,elev,gr.vgrid);

		% Copy x values so every node vertical level has an x value
		x = repmat(ob.xy.x, 1, size(sz,2)); 

		% Plot salt and sediment 
		% 2 sub plots for each variable

		%screen_sz = get(0,'ScreenSize');
		%figH = figure('Position',[ 1 screen_sz(4)/2 1200 800]);
		figH = figure('visible', 'off','color','w');
		
		%colormap(jet(256)); % Velocity
		%h(1) = subplot(2,1,1);
		%sp_pos = get(h(1), 'position'); hold on;
		%contourf(x,sz,-hvel,'linecolor','none','levelstep',0.2);
		%colorbar;
		%caxis([-2 2]);
		%buf = sprintf('RMSE %f', rmse);
		%title('Velocity m/s');
		%legH = legend('Model', 'Analytical', 'Location', 'NorthEast');
		%ylabel('Depth m');
		%plot(x(:,1:3), sz(:,1:3), 'color', [0.7 0.7 0.7], 'linewidth', 0.1);
		%axis(axisTrans);
		%set(gca,'xtick',aTTicks);
		%set(gca,'xticklabel',aTLabel);

		% Var
		h(1) = subplot(1,1,1);
		sp_pos = get(h(1), 'position'); hold on;
		contourf(x,sz,var,'linecolor','none');
		colorbar;
		caxis(axisDims); % Hardcoded with values above
		%plot(x(:,1:3), sz(:,1:3), 'color', [0.7 0.7 0.7], 'linewidth', 0.1);
		%title('Salt psu');
		%legH = legend('Model', 'Analytical', 'Location', 'NorthEast');
		axis(axisTrans);
		set(gca,'xtick',aTTicks);
		set(gca,'xticklabel',aTLabel);

		% TKE
		%h(4) = subplot(3,1,4);
		%sp_pos = get(h(4), 'position'); hold on;
		%contourf(x,sz,kine,'linecolor','none');
		%axis(axisTrans);
		%plot(x(:,1:3), sz(:,1:3), 'color', [0.7 0.7 0.7], 'linewidth', 0.1);
		%title('TKE m^2/s^2');
		%ylabel('Depth m');
		%%colormap(jet(1024));
		%%colorbar_log([0.0001 0.1]);
		%colorbar;
		%caxis([1e-6 1e-3]);
		%set(gca,'xtick',aTTicks);
		%set(gca,'xticklabel',aTLabel);

%       % Plot levels for reference
%       %plot(x(:,end), sz(:,end), 'b', 'linewidth', 0.1, 'linesmoothing', 'on'); hold on;
%       %plot(x, sz, 'color', 'w','linewidth', 0.1, 'linesmoothing', 'on'); hold on;
%       %axis([0 100000 -15 3]);

%       % Eddy diffusivity 
%       h(5) = subplot(3,1,5);
%       sp_pos = get(h(5), 'position'); hold on;
%       contourf(x,sz,tdff,'linecolor','none');
%       title('Diffusivity m^2/s');
%       axis([0 100000 -15 3]);
%       plot(x(:,1:3), sz(:,1:3), 'color', [0.7 0.7 0.7], 'linewidth', 0.1);
%       %colormap(jet(1024));
%       %colorbar_log([0.01 0.1]);
%       colorbar;
%       caxis([1e-6 1e-1]);
%       set(gca,'xtick',aTTicks);
%       set(gca,'xticklabel',aTLabel);

%       % Eddy viscosity 
%       h(6) = subplot(3,1,6);
%       sp_pos = get(h(5), 'position'); hold on;
%       contourf(x,sz,vdff,'linecolor','none');
%       title('Viscosity m^2/s');
%       axis([0 100000 -15 3]);
%       plot(x(:,1:3), sz(:,1:3), 'color', [0.7 0.7 0.7], 'linewidth', 0.1);
%       %colormap(jet(1024));
%       %colorbar_log([0.01 0.1]);
%       colorbar;
%       caxis([1e-6 1e-1]);
%       set(gca,'xtick',aTTicks);
%       set(gca,'xticklabel',aTLabel);

		% Finally add the time and text for colorbar
		% tah - text axis handle
		% ts  - timeString
		% th  - text handle
		secs = (day-1)*dt*nSteps + i*dt;
		stepTime = addtodate(startDate, secs, 'second');
		tah = axes('position', [0 0 1 1]);
		%ts  = sprintf('%2.2f %s', tideCycle, 'M_2 tidal cycles completed');
		th  = text(0.5, 0.98, datestr(stepTime, 'mm/dd/yy HH:MM:SS'), 'FontSize', 12);
		set(th, 'HorizontalAlignment', 'center');
		th  = text(0.5, 0.06, 'Distance from ocean boundary km', 'FontSize', 12);
		set(th, 'HorizontalAlignment', 'center');
		set(gca,'Visible', 'off');

		% Get ready for saves
		% n_channel_mmddyy_HHMMSS.png
		timeStr = datestr(stepTime, 'mmddyy_HHMMSS');	
		fn = sprintf('%s/%s-%s.png', plots, 'n_channel', timeStr); 
		iw = 1024;
		ih = 800;
		set(gcf,'PaperUnits','inches','PaperPosition',[0 0 iw/100.0 ih/100.0])
		print('-dpng', fn, '-r100');
		close(figH);
	end % Loop over i = 1:timeSteps - Time steps per day
	toc
end % Loop for day = 1:5 - The output days
