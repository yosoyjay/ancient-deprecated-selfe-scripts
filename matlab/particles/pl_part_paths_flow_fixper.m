function [hour] = compare_particle_paths( GR3_PATH,      ...
                                          BT_PTH_PATH,   ...
										  FT_PTH_PATH,   ...
										  HVEL_PATH,     ...
										  HOUR,          ...
										  PART_NUM,		 ...
										  varargin)      

%
%  Inputs:
%  GR3_PATH - Path to hgrid.gr3
%  BT_PTH_PATH - Path to backtrack particle.pth result file
%  FT_PTH_PATH - Path to forwardtrack particle.pth result file
%  FT_BP_PATH - Path to the particle.bp file used for forward tracking
%  HVEL_PATH - Path to n_hor_velel.64 
%  HOUR - The length of the runs
%  PART_NUM - The specific particle number to compare
%  varargin - Handle to figure if it was returned from a previous call to this function
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
%  jlopez - 12/20/2010
%
%  NOTE: No error checking is implemented at all.

% Load files / Open file for output
fg = hgrid2fg( GR3_PATH );

% Used for changing sc system to ll
% fg_ll = fg;
% [fg_ll.x, fg_ll.y] = convm2ll( fg_ll.x, fg_ll.y ); 
hor_vel = sz_readHeader( HVEL_PATH );
hour = HOUR;

[ bt_particle ft_particle ] = plot_extractParticlePaths( BT_PTH_PATH, FT_PTH_PATH, HOUR, PART_NUM )

% for centering the graph and setting an appropriate window.
% find x,y,z min and max
[ x_min x_max y_min y_max z_min z_max ] = plot_calculateAxisRange( BT_PTH_PATH, FT_PTH_PATH, HOUR, PART_NUM);

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
	set(hAxes, 'xlim', [x_min x_max], 'ylim', [y_min y_max], 'zlim', [z_max 0], ...
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
end_marker = plot3(hAxes, ft_particle(end,1), ft_particle(end,2), ft_particle(end,3),  ...
      'color', 'black',  ...
	  'Marker', 'o',     ...
	  'MarkerSize', 12);  hold on;

% Plot the paths of bt and ft run
bt_plot = line(bt_particle(:,1), bt_particle(:,2), bt_particle(:,3),  ...
      'color', 'blue',   ...
	  'LineWidth', 1,    ...
	  'Marker', '<',     ...
	  'MarkerSize', 7,   ...
	  'MarkerFaceColor', 'blue');  hold on;
ft_plot = line(ft_particle(:,1), ft_particle(:,2), ft_particle(:,3),  ...
      'color', 'black',  ...
	  'LineWidth', 1,    ...
	  'LineStyle', '--', ...
	  'Marker', '>',     ...
	  'MarkerSize', 7,   ...
	  'MarkerFaceColor', 'black'); hold on;

% Calculate distance between start of bt and end of ft run and draw it on the map
dist = calc_3d_dist( bt_particle(1,1), ft_particle(end,1), ...
                     bt_particle(1,2), ft_particle(end,2), ...
					 bt_particle(1,3), ft_particle(end,3) )

%diff_plot = line( [bt_particle(1,1) ft_particle(end,1)],  ...
%                  [bt_particle(1,2) ft_particle(end,2)],  ...
%	              [bt_particle(1,3) ft_particle(end,3)],  ...
%	              'color', 'red',                         ... 
%	              'LineWidth', 1.1);

% Annotations - Title, Legend, Hour, Difference in position
% Legend
plots_for_legend = [ start_marker end_marker bt_plot ft_plot ]; 
legend( plots_for_legend, 'Start', 'End', 'BT', 'FT', 'Location', 'NorthEastOutside');

% Title
hour = sprintf('%d', hour);
part_num = sprintf('%d', PART_NUM);
title_text = ['\bf Particle ' part_num ' after ' hour ' time step(s)' ];
title(title_text);

% Annotation - Differences in position
dist_string = sprintf('%5.2f', dist);
text_for_graph = [ 'Tracking Error: ' dist_string ' m'];
annotation('textbox', get(gca,'Position'), 'String', text_for_graph, 'FontWeight', 'bold', ...
		   'FontName', 'century gothic');

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
image_file_name = [ 'bt_ft_comp_' part_num '_' hour '.png' ];
saveas(gcf, image_file_name);
%print(gcf, image_file_name);

% Write data about distance to file
plot_data = [ part_num ' ' hour ' ' dist_string ];
fid = fopen('run_info.dat', 'a');
fprintf(fid, '%s\n', plot_data);
fclose(fid);

