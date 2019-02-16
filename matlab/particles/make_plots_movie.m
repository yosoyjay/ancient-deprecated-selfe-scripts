% This script basically just loops over a bunch of make_particle_path_plots 
% shoves them in a avi and writes the positions of each particle to a text file.
%
% Input: None
% Ouput: movie of a particles path and a text file with the position data
% Dependencies:
% compare_particle_paths.m
%
% jlopez - 12/28/2010
%
function [ ] = make_plots_movie( part_num )

videoName = [ 'plotVideo_' int2str(part_num) '.avi' ];
plotVideo = VideoWriter( videoName );
plotVideo.FrameRate = 1;
open(plotVideo);

for d=1:24
	s = ['pl_part_paths_flow_fixper(''hgrid.gr3'', ''/Users/jesse/Documents/selfe/practice/rk_6h/particle_b_1.pth'',' ...
	'''/Users/jesse/Documents/selfe/practice/rk_6h/particle_f_' int2str(d) '.pth'',' 							    ...
    '''1_hvel.64'' ,' int2str(d) ',' int2str(part_num) ')']
	plot = eval(s)
	curFrame = getframe(gcf);
	writeVideo(plotVideo, curFrame);
	close
end

close('all');


