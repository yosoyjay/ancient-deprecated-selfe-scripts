%  This script plots the paths of two runs of particle tracking to qualitatively
%  and quantitatively assess the accuracy of the tracking scheme.
%
%  Inputs:
%  GR3_PATH - Path to hgrid.gr3
%  BT_PTH_PATH - Path to backtrack particle.pth result file
%  FT_PTH_PATH - Path to forwardtrack particle.pth result file
%  FT_BP_PATH - Path to the particle.bp file used for forward tracking
%  HVEL_PATH - Path to n_hor_velel.64 
%  HOUR - The length of the runs
%  PART_NUM - The specific particle number to compare
%
%  Output
%  Well... nothing is really returned by the function... but...
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

function [hour] = compare_particle_paths( GR3_PATH,      ...
                                          BT_PTH_PATH,   ...
										  FT_PTH_PATH,   ...
										  FT_BP_PATH,    ...
										  HVEL_PATH,     ...
										  HOUR,          ...
										  PART_NUM)

% Load files / Open file for output
fg = hgrid2fg( GR3_PATH );
% Used for changing sc system to ll
% fg_ll = fg;
% [fg_ll.x, fg_ll.y] = convm2ll( fg_ll.x, fg_ll.y ); 
hor_vel = sz_readHeader( HVEL_PATH );
particles_bt = read_pth( BT_PTH_PATH );
particles_ft = read_pth( FT_PTH_PATH );
bp = read_bp ( FT_BP_PATH );
hour = HOUR;

% extract the data of a particular particle from .pth and place in an array
index = 1;

% 97-hour*4 -> i.e. 97-1*4 = 93,94,95,96 = 4 positions -> 4 positions = 1 hour
% there are 96 steps in the particle tracking.
% need to preallocate if doing large amounts of particles.
for i = (97-hour):96;
	bt_particle(index,1) = particles_bt(i,1).x(PART_NUM,1);
	bt_particle(index,2) = particles_bt(i,1).y(PART_NUM,1);
	bt_particle(index,3) = particles_bt(i,1).z(PART_NUM,1);
	index = index + 1;
end

% Need to get starting position for particle.bp to compare to end position of
% the forward tracking.  This is quite silly, ptrack should print starting position.
% Then I ignore the last position because it overshoots the backtracking by one
% time step.
ft_particle(1,1) = bp.x(PART_NUM,1);
ft_particle(1,2) = bp.y(PART_NUM,1);
ft_particle(1,3) = bp.z(PART_NUM,1);

index = 2;
for i = 2:hour;
	ft_particle(index,1) = particles_ft(i,1).x(PART_NUM,1);
	ft_particle(index,2) = particles_ft(i,1).y(PART_NUM,1);
	ft_particle(index,3) = particles_ft(i,1).z(PART_NUM,1);
    index = index + 1;
end


% For each time step, set graph parameters, draw plot, and save the figure to 
% the a file.
hg = figure('visible', 'on');

% For centering the graph and setting an appropriate window.
% Find x,y,z min and max
x_min = [ min(bt_particle(:,1)) min(ft_particle(:,1))];
x_min = min( x_min );
x_min = x_min - x_min*.0005;

x_max = [ max(bt_particle(:,1)) max(ft_particle(:,1))];
x_max = max( x_max );
x_max = x_max + x_max*.0005;

y_min = [ min(bt_particle(:,2)) min(ft_particle(:,2))];
y_min = min( y_min );
y_min = y_min - y_min*.0001;

y_max = [ max(bt_particle(:,2)) max(ft_particle(:,2))];
y_max = max( y_max );
y_max = y_max + y_max*.0001;

z_min = [ min(bt_particle(:,3)) min(ft_particle(:,3))];
z_min = min( z_min );
z_min = z_min - z_min*.0005;

z_max = [ max(bt_particle(:,3)) max(ft_particle(:,3))];
z_max = max( z_max );
z_max = z_max + z_max*.0005;


% General crap to set up figure
set(gca, 'xlim', [x_min x_max], 'ylim', [y_min y_max], 'zlim', [z_max 0]);      
plotbnd(fg);
drawelems(fg);
xlabel('Longitude (m)');
ylabel('Latitude (m)');
hold on;

% Place markers at beginning of bt run and end of ft run
start_marker = plot3(bt_particle(1,1), bt_particle(1,2), bt_particle(1,3),  ...
      'color', 'blue',   ...
	  'Marker', 'o',     ...
	  'MarkerSize', 12);  hold on;
end_marker = plot3(ft_particle(end,1), ft_particle(end,2), ft_particle(end,3),  ...
      'color', 'black',  ...
	  'Marker', 'o',     ...
	  'MarkerSize', 12);  hold on;

% Plot the paths of bt and ft run
bt_plot = plot3(bt_particle(:,1), bt_particle(:,2), bt_particle(:,3),  ...
      'color', 'blue',   ...
	  'LineWidth', 1,    ...
	  'Marker', '<',     ...
	  'MarkerSize', 7,   ...
	  'MarkerFaceColor', 'blue');  hold on;
ft_plot = plot3(ft_particle(:,1), ft_particle(:,2), ft_particle(:,3),  ...
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

diff_plot = plot3([bt_particle(1,1) ft_particle(end,1)],  ...
                  [bt_particle(1,2) ft_particle(end,2)],  ...
	              [bt_particle(1,3) ft_particle(end,3)],  ...
	              'color', 'red',                         ... 
	              'LineWidth', 1.1);

% Annotations - Title, Legend, Hour, Difference in position
% Legend
plots_for_legend = [ start_marker end_marker bt_plot ft_plot diff_plot ]; 
legend( plots_for_legend, 'Start', 'End', 'BT', 'FT', 'Diff', 'Location', 'Best');
% Title
hour = sprintf('%d', hour);
part_num = sprintf('%d', PART_NUM);
title_text = ['\bf Particle ' part_num ' after ' hour ' hour' ];
title(title_text);
% Annotation - Differences in position
dist_string = sprintf('%5.2f', dist);
text_for_graph = [ 'Difference: ' dist_string ' m'];
annotation('textbox', get(gca,'Position'), 'String', text_for_graph);

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

