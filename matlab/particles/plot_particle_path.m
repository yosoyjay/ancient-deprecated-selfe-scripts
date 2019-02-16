function [hour] = plot_particle_path( GR3_PATH,      ...
                                      BT_PTH_PATH,   ...
									  HVEL_PATH,     ...
									  HOUR,          ...
									  PART_NUM,		 ...
									  TIME_DIR ) 
%
%  Inputs:
%  GR3_PATH - Path to hgrid.gr3
%  BT_PTH_PATH - Path to backtrack particle.pth result file
%  HVEL_PATH - Path the horizonatal velocity file path
%  HOUR - The length of the runs
%  PART_NUM - The specific particle number to compare
%  TIME_DIR - time direction --- 'b' == backwards or 'f' == forwards
%
%  Output
%  hFigure - The handle to the figure to reuse for animations 
%  The resulting figure is saved in a *.png file
%  The hour and difference is also saved in a *.dat file
%
%  Dependencies:
%  hgrid2fg.m
%  read_pth.m
%  read_bp.m
%  plotbnd.m
%  drawelems.m
%  convm2ll.m
%  bcs2ll_bp - binary
%  m_elio/sz_readHeader.m
%  m_elio/sz_readTimeStep.m
%  is_valid_struct.m
%  m_elio/map_sz2hts.m
%  cal_3d_dist.m
%
%  jlopez - 2/7/2011
%
%  NOTE: No error checking is implemented at all.

% Load files / Open file for output
fg = hgrid2fg( GR3_PATH );

% Used for changing sc system to ll
fg_ll = fg;
[fg_ll.x, fg_ll.y] = convm2ll( fg_ll.x, fg_ll.y ); 
hor_vel = sz_readHeader( HVEL_PATH );
hour = HOUR;

[ bt_particle ] = plot_extractParticlePath( BT_PTH_PATH, TIME_DIR, HOUR, PART_NUM )

% for centering the graph and setting an appropriate window.
% find x,y,z min and max
[ x_min x_max y_min y_max z_min z_max ] = plot_calculateAxisRange_single( BT_PTH_PATH, TIME_DIR, HOUR, PART_NUM);

% Calculate a nice way to format the axis so they remain static for all plots
x_tick_range = abs(x_max - x_min) / 5;
y_tick_range = abs(y_max - y_min) / 5;
for i=1:6
	x_ticks(i) = x_min + x_tick_range*i;
	%x_ticks_label(i) = int2str(x_ticks(i));
	y_ticks(i) = y_min + y_tick_range*i;
	%y_ticks_label(i) = int2str(y_ticks(i));
end

if (nargin == 1)
	print 
	hFigure = varargin{1}(1);
else
	% for each time step, set graph parameters, draw plot, and save the figure to 
	% the a file.
	hFigure = figure('visible', 'on');

	% General crap to set up figure
	hAxes = axes('Parent', hFigure);
	set(hAxes, 'xlim', [x_min x_max], 'ylim', [y_min y_max], 'zlim', [z_min 0], ...
		'XLimMode', 'manual', 'YLimMode', 'manual');
	xlim = [ x_min x_max ];
	ylim = [ y_min y_max ];
	plotbnd(fg);
	drawelems(fg);
	xlabel('Longitude (m)', 'FontWeight', 'bold', 'FontName', 'century gothic');
	ylabel('Latitude (m)', 'FontWeight', 'bold', 'FontName', 'century gothic');
	hold on;
end

% Place markers at beginning of bt run and end of ft run
start_marker = plot3(hAxes, bt_particle(1,1), bt_particle(1,2), bt_particle(1,3),  ...
      'color', 'blue',   ...
	  'Marker', 'o',     ...
	  'MarkerSize', 12);  hold on;
end_marker = plot3(hAxes, bt_particle(end,1), bt_particle(end,2), bt_particle(end,3),  ...
      'color', 'black',  ...
	  'Marker', 'o',     ...
	  'MarkerSize', 12);  hold on;

% Plot the path of the run
bt_plot = line(bt_particle(:,1), bt_particle(:,2), bt_particle(:,3),  ...
      'color', 'blue',   ...
	  'LineWidth', 1,    ...
	  'Marker', '<',     ...
	  'MarkerSize', 7,   ...
	  'MarkerFaceColor', 'blue');  hold on;

%diff_plot = line( [bt_particle(1,1) ft_particle(end,1)],  ...
%                  [bt_particle(1,2) ft_particle(end,2)],  ...
%	              [bt_particle(1,3) ft_particle(end,3)],  ...
%	              'color', 'red',                         ... 
%	              'LineWidth', 1.1);

% Annotations - Title, Legend, Hour, Difference in position
% Legend
plots_for_legend = [ start_marker bt_plot ]; 
legend( plots_for_legend, 'Start', 'Path', 'Location', 'NorthEastOutside');

% Title
hour = sprintf('%d', hour);
part_num = sprintf('%d', PART_NUM);
title_text = ['\bf Particle ' part_num ' after ' hour ' time step(s)' ];
title(title_text);

% Fix up text
text_obj = findobj('type', 'text');
set(text_obj, 'fontunits', 'points');
set(text_obj, 'fontsize', 15);
set(text_obj, 'fontweight', 'bold');
set(text_obj, 'fontname', 'century gothic');

% Add velocity vectors
% Set-up stuff required for displaying flow field - from Paul's code
	it = 1;
	[d ts] = sz_readTimeStep(hor_vel, it);
	u = map_sz2hts(hor_vel, d(:,1), 1);
	v = map_sz2hts(hor_vel, d(:,2), 1);
	u = u(:,end);
	v = v(:,end);

	plot_area = [ x_min x_max y_min y_max ];

	dlat_x = (x_max - x_min) / 20;
	dlat_y = (y_max - y_min) / 20;

	[LON0, LAT0] = meshgrid([plot_area(1):dlat_x:plot_area(2)],[plot_area(3):dlat_y:plot_area(4)]);

	[LON, LAT] = meshgrid([plot_area(1):dlat_x:plot_area(2)],[plot_area(3):dlat_y:plot_area(4)]);

	utop = griddata(fg.x, fg.y, u, LON0, LAT0, 'cubic');
	vtop = griddata(fg.x, fg.y, v, LON0, LAT0, 'cubic');
	utop0 = utop; vtop0 = vtop;

	uu = interp2(LON0, LAT0, utop0, LON, LAT);
	vv = interp2(LON0, LAT0, vtop0, LON, LAT);
	% End flow field set up
	ufact = 40;
	quiver(LON, LAT, ufact*uu, ufact*vv, 0, 'k');

	it = 2;
	[d ts] = sz_readTimeStep(hor_vel, it);
	u = map_sz2hts(hor_vel, d(:,1), 1);
	v = map_sz2hts(hor_vel, d(:,2), 1);
	u = u(:,end);
	v = v(:,end);

	utop = griddata(fg.x, fg.y, u, LON0, LAT0, 'cubic');
	vtop = griddata(fg.x, fg.y, v, LON0, LAT0, 'cubic');
	utop0 = utop; vtop0 = vtop;

	uu = interp2(LON0, LAT0, utop0, LON, LAT);
	vv = interp2(LON0, LAT0, vtop0, LON, LAT);
	% End flow field set up
	ufact = 40;
	quiver(LON, LAT, ufact*uu, ufact*vv, 0, 'k', 'color', 'red');



% Save the figure - Changed to use print to enable use in script.
% 
image_file_name = [ 'pt_path_' part_num '_' hour '.png' ];
saveas(gcf, image_file_name);
%print(gcf, image_file_name);

% Write data about distance to file
fclose('all');

