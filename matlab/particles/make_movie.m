% Test to see if I can create a movie via script from the 
% cli.
%

% Load/setup grid
fg = hgrid2fg('hgrid.gr3');
fgll = fg;
[fgll.x, fgll.y] = convm2ll(fg.x, fg.y);

% Load particles
results = read_pth('particle.pth');
for i = 1,length(results)
	[results(i,1).x, results(i,1).y] = convm2ll(results(i,1).x, results(i,1).y);  
end

% Setup movie related stuf
movieObject = VideoWriter('run.avi');
open(movieObject);

% Define limits of domain for plotting
xLimits = [-124.4 -123.5];
yLimits = [46.1317 46.3186];

for i = 1:size(results,1);
	% Keep image drawing on external monitor
	set(gcf,'position', [0 0 720 400], 'color', 'w');

	% Handle elements
	plotbnd(fgll); hold on;
	%elements = drawelems(fgll); hold on;
	%dryElements = patch(fgll.x(fgll.e(1,:)),fgll.y(fgll.e(1,:)),[-.1 -.1 -.1],[.5 .5 .5]);
	%set(dryElements, 'EdgeColor', 'none');

	% Draw particles
	%scatter3(results(i,1).x, results(i,1).y, results(i,1).z, 12, 'Cdata', results(i,1).z);
	plot3(results(i,1).x, results(i,1).y, results(i,1).z, 'o');
	set(gca, 'xlim', xLimits, 'ylim', yLimits, 'zlim', [-50, 0]);

	% Set up view and capture frame
	view(0,90);
	cur_frame = getframe;
	writeVideo(movieObject, cur_frame);
end

close(movieObject);

