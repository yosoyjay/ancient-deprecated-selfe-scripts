%  This script makes an avi movie from a particle.pth file of a particular
%  paticles path.
%
%  I use it to qualitatively assess particle tracking.
%
%	TODO:
%	1. Add bathymetry - Ask Joe about how he selects/adds the shelf

% HARD CODED CRAP
GR3_PATH = 'hgrid.gr3';
BT_PTH_PATH = 'particle.pth';
FT_PTH_PATH = '../7-14-2010_1_f_24/particle.pth';
AVI_PATH = 'run.avi';
START_DATE = '07-14-2010';

% Load files / Open file for output
fg = hgrid2fg( GR3_PATH );
movie_obj = avifile( AVI_PATH );
particles_bt = read_pth( BT_PTH_PATH );
particles_ft = read_pth( FT_PTH_PATH );

% Extract the data of a particular particle from .pth and place in an array
for i = 1:96
	bt_particle(i,1) = particles_bt(i,1).x(1,1);
	bt_particle(i,2) = particles_bt(i,1).y(1,1);
	bt_particle(i,3) = particles_bt(i,1).z(1,1);
end

for i = 1:96
	ft_particle(i,1) = particles_ft(i,1).x(1,1);
	ft_particle(i,2) = particles_ft(i,1).y(1,1);
	ft_particle(i,3) = particles_ft(i,1).z(1,1);
end



% For each time step, set graph parameters, draw plot, and save the figure to 
% the .avi file.
for i = 1:size(bt_particle,1);
	hg = figure('visible', 'off');
	set(gca, 'xlim', [3.0e5 3.85e5], 'ylim', [2.5e5 3.3e5], 'zlim', [-50, 0]);      
	plotbnd(fg);
	view(0,90);
	colorbar;
	xlabel('Longitude (SPCS)');
	ylabel('Latitude (SPCS)');
	start_date = datenum( START_DATE );
	start_date = start_date - i/96;			% 96 = time steps per day (900 secs)
	annotation('textbox', get(gca,'Position'), 'String', datestr(start_date));
	hold on;
	% scatter3(results(i,1).x, results(i,1).y, results(i,1).z, 5, 'Cdata', results(20,1).z);
	comet3(bt_particle(i,1), bt_particle(i,2), bt_particle(i,3)) %'Cdata', bt_particle(i,3));
	hold on;
	movie_obj = addframe(movie_obj, hg);
end

for i = 1:size(ft_particle,1);
	comet3(ft_particle(i,1), ft_particle(i,2), ft_particle(i,3)) %'Cdata', ft_particle(i,3));
	hold on;
	movie_obj = addframe(movie_obj, hg);
end

movie_obj = close(movie_obj);
close all;

