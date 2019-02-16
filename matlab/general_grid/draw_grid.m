
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

function [hour] = draw_grid( GR3_PATH )

% Load files / Open file for output
fg = hgrid2fg( GR3_PATH );

% Used for changing sc system to ll
% fg_ll = fg;
% [fg_ll.x, fg_ll.y] = convm2ll( fg_ll.x, fg_ll.y ); 

	% General crap to set up figure
	%hAxes = axes('Parent', hFigure);
	%set(hAxes, 'xlim', [x_min x_max], 'ylim', [y_min y_max], 'zlim', [z_max 0], ...
	%	'XLimMode', 'manual', 'YLimMode', 'manual');
	plotbnd(fg);
	drawelems(fg);
	xlabel('Longitude (m)');
	ylabel('Latitude (m)');
	hold on;
end

% Set-up stuff required for displaying flow field - from Paul's code
%	it = 1;
%	[d ts] = sz_readTimeStep(hor_vel, it);
%	u = map_sz2hts(hor_vel, d(:,1), 1);
%	v = map_sz2hts(hor_vel, d(:,2), 1);
%	u = u(:,end);
%	v = v(:,end);
%
%	plot_area = [ x_min x_max y_min y_max ];
%
%	dlat_x = (x_max - x_min) / 20;
%	dlat_y = (y_max - y_min) / 20;
%
%	[LON0, LAT0] = meshgrid([plot_area(1):dlat_x:plot_area(2)],[plot_area(3):dlat_y:plot_area(4)]);
%
%	[LON, LAT] = meshgrid([plot_area(1):dlat_x:plot_area(2)],[plot_area(3):dlat_y:plot_area(4)]);
%
%	utop = griddata(fg.x, fg.y, u, LON0, LAT0, 'cubic');
%	vtop = griddata(fg.x, fg.y, v, LON0, LAT0, 'cubic');
%	utop0 = utop; vtop0 = vtop;
%
%	uu = interp2(LON0, LAT0, utop0, LON, LAT);
%	vv = interp2(LON0, LAT0, vtop0, LON, LAT);
%	% End flow field set up
%	ufact = 40;
%	% quiver(LON, LAT, ufact*uu, ufact*vv, 0, 'k');
%
%	it = 2;
%	[d ts] = sz_readTimeStep(hor_vel, it);
%	u = map_sz2hts(hor_vel, d(:,1), 1);
%	v = map_sz2hts(hor_vel, d(:,2), 1);
%	u = u(:,end);
%	v = v(:,end);
%
%	utop = griddata(fg.x, fg.y, u, LON0, LAT0, 'cubic');
%	vtop = griddata(fg.x, fg.y, v, LON0, LAT0, 'cubic');
%	utop0 = utop; vtop0 = vtop;
%
%	uu = interp2(LON0, LAT0, utop0, LON, LAT);
%	vv = interp2(LON0, LAT0, vtop0, LON, LAT);
%	% End flow field set up
%	ufact = 40;
%	% quiver(LON, LAT, ufact*uu, ufact*vv, 0, 'k', 'color', 'red');
%
%
%
%% Save the figure - Changed to use print to enable use in script.
%% 
%image_file_name = [ 'bt_ft_comp_' part_num '_' hour '.png' ];
%saveas(gcf, image_file_name);
%print(gcf, image_file_name);

% Write data about distance to file
%plot_data = [ part_num ' ' hour ' ' dist_string ];
%fid = fopen('run_info.dat', 'a');
%fprintf(fid, '%s\n', plot_data);
%fclose(fid);

