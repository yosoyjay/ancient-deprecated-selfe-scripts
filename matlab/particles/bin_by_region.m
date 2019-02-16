% This function counts the particles in each region for each time for a 
% particle.pth file and returns the results as n x 2 array.
%
% TODO: Remove hard coded region when determining particle positions
% Depends on:
% hgrid2fg.m - Reads in *.gr3 grid
% read_pth.m - Reads in *.pth particle file
% define_regions.m - Defines regions in the grid [based on db22]

function parts_by_region = bin_by_region(path_to_grid, path_to_pth)

% Obviously loads the grid and file with particles path
fg = hgrid2fg(path_to_grid);
parts = read_pth(path_to_pth);

% Get definition of regions
regions = define_regions_plus(fg);

% Each add number of particles to structure
region_names = fieldnames(regions);

% Find the number of particles in each region - HARD CODED
for i = 1:size(parts,1)
	time(i) = parts(i,1).time;
	bb_parts_index(i) = size(find(inpolygon(parts(i,1).x, parts(i,1).y, regions.baker_bay_x, regions.baker_bay_y)==1),1);
	gb_parts_index(i) = size(find(inpolygon(parts(i,1).x, parts(i,1).y, regions.grays_bay_x, regions.grays_bay_y)==1),1);
	cb_parts_index(i) = size(find(inpolygon(parts(i,1).x, parts(i,1).y, regions.calath_bay_x, regions.calath_bay_y)==1),1);
	yb_parts_index(i) = size(find(inpolygon(parts(i,1).x, parts(i,1).y, regions.youngs_bay_x, regions.youngs_bay_y)==1),1);
	em_parts_index(i) = size(find(inpolygon(parts(i,1).x, parts(i,1).y, regions.e_mouth_x, regions.e_mouth_y)==1),1);
	nc_parts_index(i) = size(find(inpolygon(parts(i,1).x, parts(i,1).y, regions.n_channel_x, regions.n_channel_y)==1),1);
	sc_parts_index(i) = size(find(inpolygon(parts(i,1).x, parts(i,1).y, regions.s_channel_x, regions.s_channel_y)==1),1);
	os_parts_index(i) = size(find(inpolygon(parts(i,1).x, parts(i,1).y, regions.or_coast_x, regions.or_coast_y)==1),1);
	ws_parts_index(i) = size(find(inpolygon(parts(i,1).x, parts(i,1).y, regions.wa_coast_x, regions.wa_coast_y)==1),1);
	us_cr_index(i) = size(find(inpolygon(parts(i,1).x, parts(i,1).y, regions.us_cr_x, regions.us_cr_y)==1),1);
end

parts_by_region = struct('time', time, 'bb', bb_parts_index, 'gb', gb_parts_index, 'cb', cb_parts_index, 'yb', yb_parts_index, 'em', ...
					  	 em_parts_index, 'nc', nc_parts_index, 'sc', sc_parts_index, 'os', os_parts_index, 'ws', ws_parts_index,        ...
					     'us', us_cr_index);


%  Testing for counting of bins
%
%sum = size(bb_parts_index) + size(gb_parts_index) + size(cb_parts_index) + ...
%	  size(yb_parts_index) + size(em_parts_index) + size(nc_parts_index) + ...
%	  size(sc_parts_index) + size(os_parts_index) + size(ws_parts_index)
%
%missing = size(parts.x) - sum
%
%region_vector = [size(bb_parts_index, 1)  size(gb_parts_index, 1)  			   ... 
%	  size(cb_parts_index, 1)   ...
%	  size(yb_parts_index, 1)   size(em_parts_index, 1)   size(nc_parts_index, 1)   ...
%	  size(sc_parts_index, 1)   size(os_parts_index, 1)   size(ws_parts_index, 1)   ...
%	  missing(1) ] 
%
