function [] = z0_comparison_animation(firstPath, secondPath) 
%  This function creates an animation comparing the sea surface
%  height from SELFE based on two different runs with different 
%  Z0 at the boundary.
%
%  It can easily be hacked up for other, similar purposes.
%
%  Input:
% 	firstPath  - Path to first run directory for modext 
%	secondPath - Path to second run directory for modext 
%
% lopezj 11/11/11

% Get data using m-elio via modext
inPath = exist('/usr/local/cmop/modext');
if inPath~=2
	addpath '/usr/local/cmop/modext'; 
end

% Dates are known a priori and hardcoded
% modext is choking on the first time step of 00:15:00s so I moved it to 00:30:00
% and decreased the loop by one
date(1) = datenum('01-01-2010 00:30:00');
for i=2:95
	date(i) = addtodate(date(i-1), 15, 'minute');
end

% Check that paths at least exist 
inPath = exist(firstPath);
if inPath~=7
	disp('firstPath not a directory.  Check for corrrectness.');	
end
inPath = exist(secondPath);
if inPath~=7
	disp('secondPath not a directory.  Check for corrrectness.');	
end

% Get data for entire domain for every time 
% Existance check only during troubleshooting.
inPath = exist('ssh_0');
if inPath~=1
	ssh_0  = modext(date, 'All', 'All', 'Surf', 'Elev', firstPath);
end
inPath = exist('ssh_03');
if inPath~=1
	ssh_03 = modext(date, 'All', 'All', 'Surf', 'Elev', firstPath);
end

% Load grid data
gridPath = sprintf('%s/hgrid.ll', firstPath);
hgrid = gr_readHGrid(gridPath);


% Create animation over every time step of extracted data 
avi = avifile('ssh_compare.avi', 'fps', 1);
for i=1:size(ssh_0.data,2)
%for i=1:2
	fig = figure('visible', 'off');
	% Create two sub-plots
	% h - handle to plots
	% ph - plot handles
	% z - SSH values to plot, this way to just change data in animation.
	%   - Grr... refreshdata does not work for patch class plots.
	h(1)  = subplot(1,2,1);
	sp_pos = get(h(1),'position');
	set(h(1), 'position', [sp_pos(1)-0.05 sp_pos(2) sp_pos(3) sp_pos(4)]);
	z1  = ssh_0.data(1,i).val;
	ph(1) = trisurf(hgrid.elem(:,3:5), hgrid.x, hgrid.y, z1, 'EdgeColor', 'none');
	set(gca,'fontsize',8);
%	set(ph(1), 'ZDataSource', 'z');
	view(0,90);
	title('Z0 = 0');
	ylim([39 50]);
	xlim([-128 -122]);

	h(2)  = subplot(1,2,2);
	sp_pos = get(h(2),'position');
	set(h(2), 'position', [sp_pos(1)+0.025 sp_pos(2) sp_pos(3) sp_pos(4)]);
	z2  = ssh_03.data(1,i).val; 
	ph(2) = trisurf(hgrid.elem(:,3:5), hgrid.x, hgrid.y, z2, 'EdgeColor', 'none');
	set(gca,'YAxisLocation', 'right');	
	set(gca,'fontsize',8);
%	set(ph(2), 'ZDataSource', 'z');
	view(0,90);
	title('Z0 = 0.3');
	ylim([39 50]);
	xlim([-128 -122]);

	% Colorbar action - Adjust position of subplots
	% hcb - handle to colorbar
	hcb = colorbar('horiz');
	set(hcb, 'Position', [0.095 0.05 0.82 0.025]);
	for j=1:2
		pos = get(h(j), 'position');
		set(h(j), 'position', [pos(1) pos(2)*1.4 pos(3)*0.98 pos(4)*0.95]);
	end
	caxis([-4 4]);
%	set(gca,'fontsize',8);
%	set(get(hcb,'ylabel'),'string','Elevation m');
	

	% Finally add the time
	% tah - text axis handle
	% ts  - timeString
	% th  - text handle
	tah = axes('position', [0 0 1 1]);
	ts  = sprintf('%s', datestr(date(i), 'mmm dd yyyy HH:MM'));
	th  = text(0.5, 0.98, ts, 'FontSize', 10);
	set(th, 'HorizontalAlignment', 'center');
	th  = text(0.04, 0.06, 'Elev. m', 'FontSize', 10);
	set(th, 'HorizontalAlignment', 'center');
	set(gca,'Visible', 'off');

%	frame = getframe(fig);
	avi = addframe(avi,fig);
	close(fig);
end
avi = close(avi); 
		
