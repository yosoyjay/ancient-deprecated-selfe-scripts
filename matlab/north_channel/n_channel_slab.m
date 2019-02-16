%r display a slab at a given depth - for R2010b matlab only
% Use in run directory. Uses files in run/outputs
% Add m-elio to the matlab path
path(path,'/usr/local/cmop/matlab/cmop/m-elio');

% Link to data files or whateverl
% ln -s /home/workspace/local0/forecasts/f22/today/run/1_salt.63 .
% ln -s /home/workspace/local0/forecasts/f22/today/run/1_zcor.63 .
% ln -s /home/workspace/local0/forecasts/f22/today/run/hgrid.gr3 .

% Read a lat/lon .ll grid for other coordinates
gr.hgrid=gr_readHGrid('hgrid.gr3');

movie = avifile('bl.avi')
fig = figure('visible', 'off', 'color', 'w'); 

% Time steps
startTime = datenum('04/23/2012 00:15:00')
startStep = 1;
skipStep = 2;
endStep = 96;
startDay = 1;
endDay = 25;
times = startTime:((skipStep*15)/86400):startTime+endDay-1;

for day = startDay:endDay
	% Read the header for variable and vertical grid
	f = sprintf('outputs/%d_salt.63', day)
	h = sz_readHeader(f);
	f = sprintf('outputs/%d_zcor.63', day)
	hz = sz_readHeader(f);

	for it = startStep:skipStep:endStep
		[d ts]=sz_readTimeStep(h,it); 
		[dz tsz]=sz_readTimeStep(hz,it); 

		% Variable of interest
		[s] = map_sz2hts(h,d(:),1);
		% Vertical grid as a .63 file
		[z] = map_sz2hts(hz,dz(:),1);

		% also can do this:
		[s0] = s(:,end); % End is the surface layer
		% or this:
		vbot = 18; % Based on vgrid.  Last z-level, first s-level
		[sn] = s(:,vbot); % vbot is the bottom layer for sigma layers

		% Extract the slab
		%st=filter_depth_q(s, z, myslab);

		% set limits
% North Channel Focus
%		xlim([335120-200 356070+200]);
%		ylim([289380-200 295350+200]);
% Estuary
%		xlim([330000 360000]);
%		ylim([280000 300000]);
		c = [0 32];

		% Plot surface
		subplot(2,1,1);
		gr_plot(gr.hgrid,s0,c);
		xlim([330000 360000]);
		ylim([280000 300000]);
    ylabel('Latitude (m)');
    xlabel('Longitude (m)');
    set(gca,'FontSize', 14);


		% Plot bottom
		subplot(2,1,2);
		xlim([330000 360000]);
		ylim([280000 300000]);
	  gr_plot(gr.hgrid,sn,c);
	  ylabel('Latitude (m)');
    xlabel('Longitude (m)');
    set(gca,'FontSize', 14);


		%title
		%idx = (day-1)*endStep + it;
    %title(datestr(times(idx), 'mm/dd/yyyy HH:MM:SS'));
		%set(gca,'FontSize', 14);

		movie = addframe(movie, fig);
		clf(fig);
	end
end
movie = close(movie)
