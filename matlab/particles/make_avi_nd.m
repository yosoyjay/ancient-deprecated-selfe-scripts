%  This script makes an avi movie from a particle.pth file.
%
%	TODO:
%	1. Add bathymetry 

% HARD CODED CRAP
GR3_PATH = 'hgrid.gr3';
PTH_PATH = 'particle.pth';
AVI_PATH = 'run.avi';
START_DATE = '07-07-2010';


% Load files / Open file for output
fg = hgrid2fg( GR3_PATH );
results = read_pth( PTH_PATH );
movie_obj = avifile( AVI_PATH );

% For each time step, set graph parameters, draw plot, and save the figure to 
% the .avi file.
for i = 1:size(results,1);
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
	scatter3(results(i,1).x, results(i,1).y, results(i,1).z, 10, 'Cdata', results(20,1).z);
	movie_obj = addframe(movie_obj, hg);
end

movie_obj = close(movie_obj);
close all;

