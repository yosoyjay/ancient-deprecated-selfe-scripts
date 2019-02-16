function [] = make_plots_anim(grHGrid,data,anim,startDate,type) 
% [] = make_plots_anim(grHGrid,data,anim,startDate) 
%
% Creates a set of plots with option of animation of 2d data of model output.
%
% Input:
%	grHGrid - HGrid in style used with m-elio library.
%	data	- Data extract from model, I assume it to have the following shape
%			  (day, timeStep, nodeData).
%	anim	- Create animation? 1 = yes, 0 = no.
%	startDate - Start date for identifing time of plot (Date num or format that can
%				easily be converted to datenum.
%	type - 1 = instantaneous, 2 = hour, 3 = region for kflow
%
% Output:
%	plots   - Plots for every time step
%	anim	- If selected an animation of the plots
%
%

% M-elio library used
addpath '/usr/local/cmop/matlab/cmop/m-elio/';

% Constants used in plotting
skipStep = 1  % Number of time steps to skip in dim 2 of the data
startDay = 3  % Index in dim 1 of data to start with
endDay = 10    % Index in dim 1 of data to end with
figPos = [1 1 1200 800];
fontSize = 14;

xlim = [3.8e+5 5e+5];		% River only
ylim = [1.9e+5 2.9e+5];		% River only
if type == 1
	cbar = [0 0.2];
elseif type == 2
	cbar = [0 40];
else
 	cbar = [0 20e+9];				% K_flux [0 0.2] instantaneous [0 800] hour [0 20*e9] region
end
% Corresponding colorbar argument commented out on line 52;

% Set up movie if required
if anim == 1
	hAnim = avifile('plot_anim.avi');
end

% Datenum for dates?
date = datenum(startDate);

% Loop over days
%for day = 1:size(data,1)
for day = startDay:endDay
	% Loop over time steps
	for ts = 1:skipStep:size(data,2)
		% Create figure
		hFig = figure('visible','off');
		subplot(1,1,1);
		set(hFig, 'color','w');

		% Need to reshape data to Mx1 matrix and assume colobar range [0 0.2]
		gr_plot(grHGrid, reshape(data(day,ts,:),size(data,3),1), cbar);	
		%gr_plot(grHGrid, reshape(data(day,ts,:),size(data,3),1));	

		% Annotate and zoom in and such
		set(hFig, 'Position', figPos)
		set(gca,'xlim', xlim);
		set(gca,'ylim', ylim);
		grid on;
		box on;
	
		% Axis labels
		set(gca,'FontSize', fontSize);
		xlabel('X - SPCS (m)');
		ylabel('Y - SPCS (m)');

		% Title
		%date = addtodate(date,30,'minute');
		date = addtodate(date,1,'hour');
		str = sprintf('k_{flux} %s', datestr(date,'mmm dd, yyyy HH:MM:SS'));
		title(str);

		% Colorbar
		hCBar = colorbar;
		set(gca,'FontSize', fontSize);
		set(get(hCBar,'ylabel'), 'String', 'cm h^{-1}');

		% Get ready for saves
		fName = sprintf('k_flux_%d_%d.png', day, ts);
		width = 1024;
		height = 800;
		set(ehFig,'PaperUnits','inches','PaperPosition',[0 0 width/100.0 height/100.0]);
		print('-dpng', fName,'-r100');

		% Save frame to avi if required
		if anim == 1
			hAnim = addframe(hAnim, hFig);
		end
		close(hFig);
	end
end

if anim == 1
	hAnim = close(hAnim);
end
