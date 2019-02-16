% Rough sketch of begining data anlysis of particle tracking
% Idea:
% 1.  Load particles from backtrack and forward tracking
% 2.  Compare origin and destination
% 3.  Compare time from origin to destination
% 4.  Compare delta from backword -> forward & forward -> backword
% Load particle tracking results
% bt_particles = partpath(PATH_TO_BACKTRACKING_RESULTS);
% ft_particles = partpath(PATH_TO_FOWARDTRACKING_RESULTS);
%
%
% TODO: Remove hard coded region when determining particle positions

cd '/Users/jesse/Documents/selfe/practice/';

% Obviously loads the grid and file with particles path
fg = hgrid2fg('hgrid.gr3');
parts = read_pth('particle.pth');

% Get definition of regions
regions = define_regions(fg);

% Each add number of particles to structure
region_names = fieldnames(regions);

% Find the number of particles in each region - HARD CODED
for i = 1:size(parts,1)
	bb_parts_index(i) = size(find(inpolygon(parts(i,1).x, parts(i,1).y, regions.baker_bay_x, regions.baker_bay_y)==1),1);
	gb_parts_index(i) = size(find(inpolygon(parts(i,1).x, parts(i,1).y, regions.grays_bay_x, regions.grays_bay_y)==1),1);
	cb_parts_index(i) = size(find(inpolygon(parts(i,1).x, parts(i,1).y, regions.calath_bay_x, regions.calath_bay_y)==1),1);
	yb_parts_index(i) = size(find(inpolygon(parts(i,1).x, parts(i,1).y, regions.youngs_bay_x, regions.youngs_bay_y)==1),1);
	em_parts_index(i) = size(find(inpolygon(parts(i,1).x, parts(i,1).y, regions.e_mouth_x, regions.e_mouth_y)==1),1);
	nc_parts_index(i) = size(find(inpolygon(parts(i,1).x, parts(i,1).y, regions.n_channel_x, regions.n_channel_y)==1),1);
	sc_parts_index(i) = size(find(inpolygon(parts(i,1).x, parts(i,1).y, regions.s_channel_x, regions.s_channel_y)==1),1);
	os_parts_index(i) = size(find(inpolygon(parts(i,1).x, parts(i,1).y, regions.or_coast_x, regions.or_coast_y)==1),1);
	ws_parts_index(i) = size(find(inpolygon(parts(i,1).x, parts(i,1).y, regions.wa_coast_x, regions.wa_coast_y)==1),1);
end
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
