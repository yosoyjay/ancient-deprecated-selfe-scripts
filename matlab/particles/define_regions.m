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

x_min=3.31e5;
x_max=3.68e5;
y_min=2.71e5;
y_max=3.09e5;

g_index = find(grid_bnd_x > x_min & grid_bnd_x < x_max & ...
			 grid_bnd_y > y_min & grid_bnd_y < y_max);

%
% Define different regions in the grid.
%

% 
% Region: Baker Bay
% 
baker_bay_x = [grid_bnd_x(g_index(461:514)); grid_bnd_x(g_index(461))];
baker_bay_y = [grid_bnd_y(g_index(461:514)); grid_bnd_y(g_index(461))];

% 
% Region: Grays Bay
% 
grays_bay_x = [grid_bnd_x(g_index(361:408)); grid_bnd_x(g_index(361))];
grays_bay_y = [grid_bnd_y(g_index(361:408)); grid_bnd_y(g_index(361))];

% 
% Region: Calathmet Bay
% 
calath_bay_x = [grid_bnd_x(g_index(266:350)); grid_bnd_x(g_index(350)); ...
			    grid_bnd_x(g_index(266))];
calath_bay_y = [grid_bnd_y(g_index(266:350)); grid_bnd_y(g_index(266)); ...
			    grid_bnd_y(g_index(266))];

% 
% Region: Youngs Bay
% 
youngs_bay_x = [grid_bnd_x(g_index(119:219)); grid_bnd_x(g_index(119))]; 
youngs_bay_y = [grid_bnd_y(g_index(119:219)); grid_bnd_y(g_index(119))]; 

% 
% Region: Estuary Mouth 
% 
e_mouth_x = [grid_bnd_x(g_index(47:96)); grid_bnd_x(g_index(461)); ... 
			 grid_bnd_x(g_index(514)); grid_bnd_x(g_index(514:533)); ...
			 grid_bnd_x(g_index(47))];
e_mouth_y = [grid_bnd_y(g_index(47:96)); grid_bnd_y(g_index(461)); ...
			 grid_bnd_y(g_index(514)); grid_bnd_y(g_index(514:533)); ...
			 grid_bnd_y(g_index(47))];

% 
% Region: North Channel 
% 
n_channel_x = [grid_bnd_x(g_index(351:361)); grid_bnd_x(g_index(408)); ... 
		grid_bnd_x(g_index(408:461)); grid_bnd_x(g_index(101)); ...
		grid_bnd_x(g_index(220)); grid_bnd_x(g_index(759:765)); ...
		grid_bnd_x(g_index(350)); grid_bnd_x(g_index(351))];

n_channel_y = [grid_bnd_y(g_index(351:361)); grid_bnd_y(g_index(408)); ...
		grid_bnd_y(g_index(408:461)); grid_bnd_y(g_index(76));  ...
		grid_bnd_y(g_index(90)); grid_bnd_y(g_index(759:765));  ...
		grid_bnd_y(g_index(765)); grid_bnd_y(g_index(351))];

% 
% region: South Channel 
% 
s_channel_x = [grid_bnd_x(g_index(96:119)); grid_bnd_x(g_index(219));  ...
		grid_bnd_x(g_index(219:266)); grid_bnd_x(g_index(350)); ...
		grid_bnd_x(g_index(350)); grid_bnd_x(g_index(765));		...
		grid_bnd_x(g_index(765:776)); grid_bnd_x(g_index(759)); ...
		grid_bnd_x(g_index(220)); ...
		grid_bnd_x(g_index(101)); grid_bnd_x(g_index(96))];

s_channel_y = [grid_bnd_y(g_index(96:119)); grid_bnd_y(g_index(219));  ...
		grid_bnd_y(g_index(219:266)); grid_bnd_y(g_index(266)); ...
		grid_bnd_y(g_index(765)); grid_bnd_y(g_index(765));		... 
		grid_bnd_y(g_index(765:776)); grid_bnd_y(g_index(759)); ... 
		grid_bnd_y(g_index(90));  ...
		grid_bnd_y(g_index(76)); grid_bnd_y(g_index(96))];    

%
% region: Upstream CR
%

% 
% region: Oregon Shelf 
%
% TODO: Figure out why grid_bnd_y(g_index(450)) works.
or_coast_x = [grid_bnd_x(g_index(1:47)); grid_bnd_x(g_index(45));  ...
		grid_bnd_y(g_index(450)); grid_bnd_y(g_index(450));  ...
		grid_bnd_x(g_index(1))];

or_coast_y = [grid_bnd_y(g_index(1:47)); grid_bnd_y(g_index(461)); ...
		grid_bnd_y(g_index(450)); grid_bnd_y(g_index(1));	 ...
		grid_bnd_y(g_index(1))];

% 
% Region: Washington Shelf 
%
% TODO: Figure out why grid_bnd_y(g_index(450)) works.
wa_coast_x = [grid_bnd_x(g_index(45)); grid_bnd_x(g_index(533:560));  ...
		grid_bnd_y(g_index(450)); grid_bnd_y(g_index(450));     ...
		grid_bnd_x(g_index(45))];

wa_coast_y = [grid_bnd_y(g_index(461)); grid_bnd_y(g_index(533:560)); ...
		grid_bnd_y(g_index(560)); grid_bnd_y(g_index(450));     ...
		grid_bnd_y(g_index(461))];

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
set(gca,'xlim',[3.1e5 3.85e5], 'ylim', [2.6e5 3.3e5]);

regions = struct( 'baker_bay_x', baker_bay_x, 'baker_bay_y', baker_bay_y, ...
				  'grays_bay_x', grays_bay_y, 'grays_bay_y', grays_bay_y, ...
				  'calath_bay_x', calath_bay_x, 'calath_bay_y', calath_bay_y, ...
				  'youngs_bay_x', youngs_bay_x, 'youngs_bay_y', youngs_bay_y, ...
				  'e_mouth_x', e_mouth_x, 'e_mouth_y', e_mouth_y, 		  ...
				  'n_channel_x', n_channel_x, 'n_channel_y', n_channel_y, ...
				  's_channel_x', s_channel_x, 's_channel_y', s_channel_y, ...
				  'or_coast_x', or_coast_x, 'or_coast_y', or_coast_y,     ...
				  'wa_coast_x', wa_coast_x, 'wa_coast_y', wa_coast_y );

