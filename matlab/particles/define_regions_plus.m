% This fun_channeltion defines the regions in the estuary and ocean used for particle 
% tracking. The fun_channeltion is (not overloaded) to accept a *.gr3 file or accept 
% the path to a *.gr3 file.  Returns a cell array with each region's name definition 
% as a structure and each region definiton as two vectors (x,y) defined in
% terms of nodes. 
%
% jlopez 9/6/2010
% TODO: Implement using this fun_channeltion with path to fg.

% fg = Your *.gr3 file
function [regions] = define_regions(fg)

% Define the region under inspection.
% Load the boundary of the grid to ease definition of regions by 
% limiting the number of nodes to deal with.  
% Then limit the boundary to the part of the grid to the maximum
% extend of salinity intrustion to just off the shelf.
grid_bnd_x = fg.x(fg.bnd(:,1));
grid_bnd_y = fg.y(fg.bnd(:,1));

x_min=3.30e5;
x_max=4.10e5;
y_min=2.6e5;
y_max=3.3e5;

g_index = find(grid_bnd_x > x_min & grid_bnd_x < x_max & ...
			 grid_bnd_y > y_min & grid_bnd_y < y_max)

%
% Define different regions in the grid.
%

% 
% Region: Baker Bay - 53, 1  
% 
baker_bay_x = [grid_bnd_x(g_index(840:893)); grid_bnd_x(g_index(840))];
baker_bay_y = [grid_bnd_y(g_index(840:893)); grid_bnd_y(g_index(840))];

% 
% Region: Grays Bay - 47, 1
% 
grays_bay_x = [grid_bnd_x(g_index(740:787)); grid_bnd_x(g_index(740))];
grays_bay_y = [grid_bnd_y(g_index(740:787)); grid_bnd_y(g_index(740))];

% 
% Region: Calathmet Bay - 34, 1, 1
% 
calath_bay_x = [grid_bnd_x(g_index(283:422)); grid_bnd_x(g_index(350)); ...
			    grid_bnd_x(g_index(283))];
calath_bay_y = [grid_bnd_y(g_index(283:422)); grid_bnd_y(g_index(283)); ...
			    grid_bnd_y(g_index(283))];

% 
% Region: Youngs Bay - 100, 1
% 
youngs_bay_x = [grid_bnd_x(g_index(135:237)); grid_bnd_x(g_index(135))]; 
youngs_bay_y = [grid_bnd_y(g_index(135:237)); grid_bnd_y(g_index(135))]; 

% 
% Region: Estuary Mouth - 49, 1, 1, 19  
% 
e_mouth_x = [grid_bnd_x(g_index(62:113)); grid_bnd_x(g_index(840)); ... 
			 grid_bnd_x(g_index(893)); grid_bnd_x(g_index(893:912)); ...
			 grid_bnd_x(g_index(62))];
e_mouth_y = [grid_bnd_y(g_index(62:113)); grid_bnd_y(g_index(840)); ...
			 grid_bnd_y(g_index(893)); grid_bnd_y(g_index(893:912)); ...
			 grid_bnd_y(g_index(62))];

% 
% Region: North Channel - 10, 1, 53, 1, 1, 6, 1, 1
% 
n_channel_x = [grid_bnd_x(g_index(702:740)); grid_bnd_x(g_index(787)); ... 
			   grid_bnd_x(g_index(787:840)); grid_bnd_x(g_index(120)); ...
			   grid_bnd_x(g_index(250)); grid_bnd_x(g_index(1344:1351)); ...
			   grid_bnd_x(g_index(1317)); grid_bnd_x(g_index(413));    ...
			   grid_bnd_x(g_index(702))];
n_channel_y = [grid_bnd_y(g_index(702:740)); grid_bnd_y(g_index(787)); ...
			   grid_bnd_y(g_index(787:840)); grid_bnd_y(g_index(91));  ...
			   grid_bnd_y(g_index(100)); grid_bnd_y(g_index(1344:1351)); ...
			   grid_bnd_y(g_index(1317)); grid_bnd_y(g_index(1325));   ...
			   grid_bnd_y(g_index(702))];

% 
% region: South Channel - 23, 1, 47, 1, 1, 1, 11, 1, 1, 1, 1
% 
s_channel_x = [grid_bnd_x(g_index(120)); grid_bnd_x(g_index(113:135));  ...
			   grid_bnd_x(g_index(237:283)); grid_bnd_x(g_index(350)); ...
			   grid_bnd_x(g_index(422));  ... 
			   grid_bnd_x(g_index(413));  ...
			   grid_bnd_x(g_index(1317)); grid_bnd_x(g_index(1351:1361)); ...
			   grid_bnd_x(g_index(1344)); grid_bnd_x(g_index(250));		...
			   grid_bnd_x(g_index(120))];

s_channel_y = [grid_bnd_y(g_index(91)); grid_bnd_y(g_index(113:135));   ...
			   grid_bnd_y(g_index(237:283)); grid_bnd_y(g_index(283));  ...
			   grid_bnd_y(g_index(422));  ... 
			   grid_bnd_y(g_index(1325)); ...
			   grid_bnd_y(g_index(1317)); grid_bnd_y(g_index(1351:1361)); ...
			   grid_bnd_y(g_index(1344)); grid_bnd_y(g_index(100));		...
			   grid_bnd_y(g_index(91))];
%
% region: Upstream CR
%
% TODO: Add end of island to domain 
us_cr_x = [grid_bnd_x(g_index(422:702)); grid_bnd_x(g_index(422))];
us_cr_y = [grid_bnd_y(g_index(422:702)); grid_bnd_y(g_index(422))]; 

% 
% region: Oregon Shelf  - 57, 1, 1, 1, 1
%
% TODO: Figure out why grid_bnd_y(g_index(450)) works.
or_coast_x = [grid_bnd_x(g_index(5:62)); grid_bnd_x(g_index(910));  ...
		grid_bnd_y(g_index(450)); grid_bnd_y(g_index(450));  ...
		grid_bnd_x(g_index(5))];


or_coast_y = [grid_bnd_y(g_index(5:62)); grid_bnd_y(g_index(820)); ...
		grid_bnd_y(g_index(820)); grid_bnd_y(g_index(5));	 ...
		grid_bnd_y(g_index(5))];

% 
% Region: Washington Shelf - 1, 27, 1, 1, 1, 
%
% TODO: Figure out why grid_bnd_y(g_index(450)) works.
wa_coast_x = [grid_bnd_x(g_index(910)); grid_bnd_x(g_index(912:960));  ...
		grid_bnd_y(g_index(450)); grid_bnd_y(g_index(450));     ...
		grid_bnd_x(g_index(910))];

wa_coast_y = [grid_bnd_y(g_index(820)); grid_bnd_y(g_index(912:960)); ...
		grid_bnd_y(g_index(960)); grid_bnd_y(g_index(820));     ...
		grid_bnd_y(g_index(820))];

% Plot the region
clf;
plotbnd(fg);
bb_line = line(baker_bay_x, baker_bay_y);
gb_line = line(grays_bay_x, grays_bay_y);
cb_line = line(calath_bay_x, calath_bay_y);
yb_line = line(youngs_bay_x, youngs_bay_y);
em_line = line(e_mouth_x, e_mouth_y);
n_channel_line = line(n_channel_x, n_channel_y);
s_channel_line = line(s_channel_x, s_channel_y);
or_coast_line = line(or_coast_x, or_coast_y);
wa_coast_line = line(wa_coast_x, wa_coast_y);
us_cr_line = line(us_cr_x, us_cr_y);
set(gca,'xlim',[3.1e5 3.85e5], 'ylim', [2.6e5 3.3e5]);

regions = struct( 'baker_bay_x', baker_bay_x, 'baker_bay_y', baker_bay_y, ...
				  'grays_bay_x', grays_bay_y, 'grays_bay_y', grays_bay_y, ...
				  'calath_bay_x', calath_bay_x, 'calath_bay_y', calath_bay_y, ...
				  'youngs_bay_x', youngs_bay_x, 'youngs_bay_y', youngs_bay_y, ...
				  'e_mouth_x', e_mouth_x, 'e_mouth_y', e_mouth_y, 		  ...
				  'n_channel_x', n_channel_x, 'n_channel_y', n_channel_y, ...
				  's_channel_x', s_channel_x, 's_channel_y', s_channel_y, ...
				  'or_coast_x', or_coast_x, 'or_coast_y', or_coast_y,     ...
				  'wa_coast_x', wa_coast_x, 'wa_coast_y', wa_coast_y,	  ...
				  'us_cr_x', us_cr_x, 'us_cr_y', us_cr_y);

