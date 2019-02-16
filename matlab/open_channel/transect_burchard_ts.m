function [] = plotTransect(runPath, startDay, endDay, xPos, levels)
% This function creates a plot for a transect to emulate those in Burchard et al. 2005 
%
% Input:     
%   runPath - Path to run directory, assumes outputs are already combined in outputs dir
%           - Assumes *.bp file named burchard.bp in directory made specifically for this.
%			- A bunch of stuff is hardcoded based on this file so be aware
%			- Time step # (Not time) to plot data? 
%	startDay - Day number of run (File number 1_hvel, 2_hvel, etc...)
%   endDay   - Day number of run (File number 6_hvel, 12_hvel, etc...)
%	xPos    - Index in burchard.bp that give the x position of where I want to make time series
%   levels  - Index of the s-levels to plot the data
%   movie   - 1 - Makes a avi, 0 - No avi
%
% Output:
%
% lopezj - 12/14/11
tic;

% Constants etc.
addpath('/usr/local/cmop/matlab/cmop/m-elio');
addpath('/home/workspace/users/lopezj/scripts/matlab');
xPos = 101; % Place in burchard.bp that I want to get time series.

% Make a plots folder if it doesn't exists for the output 
plots = sprintf('%s/plots', runPath);
if (~exist(plots,'file'))
	mkdir(plots);
end

% Build points file
bpPath = sprintf('%s/burchard.bp', runPath);
tr.hgrid = gr_readHGrid(bpPath);

% Get information about the number of days in run from param.in
% Assumes it is on line 151
paramPath = sprintf('%s/param.in', runPath);
param = fopen(paramPath, 'r');
C = textscan(param,'%s','delimiter','\n');
dayRawStr = C{1,1}(151);
dayRawStr = regexp(dayRawStr,' ','split');
days = str2num(cell2mat(dayRawStr{1,1}(3)));

% Prepare arrays for data reogranization
hvelPath = sprintf('%s/outputs/%d_hvel.64', runPath, startDay);
hHvel = sz_readHeader(hvelPath);
nSteps = hHvel.nSteps;
hvel = zeros((endDay-startDay+1)*nSteps,max(size(levels)));
salt = zeros((endDay-startDay+1)*nSteps,max(size(levels)));
kine = zeros((endDay-startDay+1)*nSteps,max(size(levels)));
tdff = zeros((endDay-startDay+1)*nSteps,max(size(levels)));
sed = zeros((endDay-startDay+1)*nSteps,max(size(levels)));

for j=startDay:endDay
	% Deal with h- and v-grid
	% Note that sz layers must be extracted from binary output file headers
	% Open up files for all variables I'm looking at in open channel
	hgridPath = sprintf('%s/hgrid.gr3', runPath);
	gr.hgrid = gr_readHGrid(hgridPath);
	hvelPath = sprintf('%s/outputs/%d_hvel.64', runPath, j);
	hHvel = sz_readHeader(hvelPath);
	kinePath = sprintf('%s/outputs/%d_kine.63', runPath, j);
	hKine = sz_readHeader(kinePath);
	tdffPath = sprintf('%s/outputs/%d_tdff.63', runPath, j);
	hTdff = sz_readHeader(tdffPath);
	vdffPath = sprintf('%s/outputs/%d_vdff.63', runPath, j);
	hVdff = sz_readHeader(vdffPath);
	elevPath = sprintf('%s/outputs/%d_elev.61', runPath, j);
	hElev = sz_readHeader(elevPath);
	saltPath = sprintf('%s/outputs/%d_salt.63', runPath, j);
	hSalt = sz_readHeader(saltPath);
	%sedPath = sprintf('%s/outputs/%d_trcr_1.63', runPath, j);
	%hSed = sz_readHeader(sedPath);
	gr.vgrid = hKine.vgrid;

	% Compute transect
	[ob]= ob_ini_fromTrasect(gr, bpPath);
	trLen = cumsum(sqrt((ob.xy.x(2:end)-ob.xy.x(1:end-1)).^2 + ...
						(ob.xy.y(2:end)-ob.xy.y(1:end-1)).^2));
	trLen = [0; trLen];
	nSteps = hHvel.nSteps;

	% Collect data from every time step
	for i=1:nSteps
		% Read timestep of data
		[hvelData varTs] = sz_readTimeStep(hHvel,i);
		[kineData varTs] = sz_readTimeStep(hKine,i);
		[tdffData varTs] = sz_readTimeStep(hTdff,i);
		[vdffData varTs] = sz_readTimeStep(hVdff,i);
		[saltData varTs] = sz_readTimeStep(hSalt,i);
		%[sedData varTs]  = sz_readTimeStep(hSed,i);

		% Map sz levels to depths
		dHvel = map_sz2hts(hHvel, hvelData(:,1));
		dKine = map_sz2hts(hKine, kineData);
		dTdff = map_sz2hts(hTdff, tdffData);
		dVdff = map_sz2hts(hVdff, vdffData);
		dSalt = map_sz2hts(hSalt, saltData);
		%dSed  = map_sz2hts(hSed, sedData);

		% Do this and you just end up with the data that you want in transect
		% ob is the object that has the transect information
		temp_hvel = ob.xy.H*double(dHvel);
		temp_kine = ob.xy.H*double(dKine);
		temp_tdff = ob.xy.H*double(dTdff);
		temp_vdff = ob.xy.H*double(dVdff);
		temp_salt = ob.xy.H*double(dSalt);
		%temp_sed  = ob.xy.H*double(dSed);

		% Construct vertical grid
		% Read timestep of elevations
		elevData = sz_readTimeStep(hElev,i);
		elev     = ob.xy.H*double(elevData);
		depths   = ob.xy.H*gr.hgrid.depth;
		sz=sz_computeZlevels(depths,elev,gr.vgrid);

		% Copy x values so every node vertical level has an x value
		x = repmat(ob.xy.x, 1, size(sz,2)); 

		hvel((j-1)*nSteps+nSteps,:) = temp_hvel(1,levels);   
		kine((j-1)*nSteps+nSteps,:) = temp_kine(1,levels);   
		salt((j-1)*nSteps+nSteps,:) = temp_salt(1,levels);   
		vdff((j-1)*nSteps+nSteps,:) = temp_vdff(1,levels);
		%sed((j-1)*nSteps+nSteps,:) = temp_sed(1,levels);   

	end % Loop over i = 1:timeSteps - Time steps per day
end	% Loop for j = 1:5 - The output days

% Axis limits for all variables
axisDimsTdff = [10e-6 0];
axisDimsVdff = [10e-6 0];
axisDimsHvel = [-2 2];
axisDimsKine = [0 0.15];
axisDimsSed  = [0.1 1.4];
axisDimsSalt = [0 30];

% Colors
color(1) = 'r';
color(2) = 'k';
color(3) = 'b';
color(4) = 'g';

% Plot data and analytical solution
% 4 sub plots for each variable
figH = figure('color','w');

% Velocity
h(1) = subplot(1,4,1);
sp_pos = get(h(1), 'position'); hold on;
for lev = 1:size(levels)
	plot(1:(nSteps*(endDay-startDay+1)),hvel(:,lev),'color',color(lev));
end
ylabel('Velocity m/s');

% Salt
h(2) = subplot(1,4,2);
sp_pos = get(h(2), 'position'); hold on;
for lev = 1:size(levels)
	plot(1:(nSteps*(endDay-startDay+1)),salt(:,lev),'color',color(lev));
end
ylabel('Salinity psu');

% Vdff 
h(3) = subplot(1,4,3);
sp_pos = get(h(3), 'position'); hold on;
for lev = 1:size(levels)
	plot(1:(nSteps*(endDay-startDay+1)),vdff(:,lev),'color',color(lev));
end
ylabel('Vdff m^2/s');

% TKE
h(4) = subplot(1,4,4);
sp_pos = get(h(4), 'position'); hold on;
for lev = 1:size(levels)
	plot(1:(nSteps*(endDay-startDay+1)),hvel(:,lev),'color',color(lev));
end
ylabel('TKE m^2/s^2');


% Finally add the time and text for colorbar
% tah - text axis handle
% ts  - timeString
% th  - text handle
secs = (j-1)*900*nSteps + i*900;
tideCycle = secs/(24.84*60*60);
tah = axes('position', [0 0 1 1]);
ts  = sprintf('%2.2f %s', tideCycle, 'M_2 tidal cycles completed');
th  = text(0.5, 0.98, ts, 'FontSize', 12);
set(th, 'HorizontalAlignment', 'center');
th  = text(0.5, 0.06, 'Distance from ocean boundary km', 'FontSize', 12);
set(th, 'HorizontalAlignment', 'center');
set(gca,'Visible', 'off');

% Get ready for saves
fn = sprintf('%s/%s-%06d.png', plots, 'burchard', secs); 
iw = 1024;
ih = 800;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 iw/100.0 ih/100.0])
print('-dpng', fn, '-r100');
