% Matlab masterscript for generating and plotting particle tracking forecasts
% (all necessary m-files are assumed to be available in the same file as the
% masterscript). In the case of positions, data may be supplied in either
% Lat/Lon (decimal degrees or deg decmin) or x/y pairs. This script was
% developed and tested to run on FCAST00.


% Set all user-defined options here...

% This is the name of the forecast as it appears in the directory call ('dev', 'db16', etc.)
modname='dev';

% The approximate number of particles desired per release. If relshape only
% includes one point, this value will default to 1. If relshape includes 2
% points (a transect), this number defines how many evenly spaced release
% locations there are along the transect. If relshape includes 3+
% positions, partnum defines how many particles will be released as a
% cluster from each release position
%(total particles = partnum*releasepositions*repetitions)
partnum=40;

% Defines the region in which particles are released. If only 1 point is
% listed, the particle will be released there; 2 points makes a release transect; 
% 3+ points will create a composite plot showing the projected trajectories for each
% point. NOTE: When creating composite plots, avoid choosing release points that
% lie withing 500 meters of the coastline.

% relshape = [46.23039 -124.10988;
%             46.23039 -124.13495; 
%             46.23039 -124.15707; 
%             46.26430 -124.10988; 
%             46.26430 -124.13495; 
%             46.26430 -124.15707];

% relshape = [46 15.3661 -124 17.2179];
% relshape = [46.245 -124.248];
% relshape = [46 16.6362 -124 15.8721];
% relshape = [46 14.0779 -124 16.8479];

% relshape = [46.23039 -124.15707;
%              46.26725 -124.15117;
%              46.23039 -124.15117;
%              46.24513 -124.12905;
%              46.26135 -124.14527];

relshape = [46.19942 -124.10545;
            46.19942 -124.13200; 
            46.19942 -124.16297; 
            46.19942 -124.19837];

% When to release particles (year, month, day, hour (PST), min . Multiple releases
% are created by adding new rows with the same format. Currently is not set up to
% span 2 forecast days. Note: reltime should only include multiple release
% times when building composite plots (i.e., relshape includes 3+ positions)
% reltime=[2008 9 25 17 45];
reltime=[2010 7 6 00 15];

%_____________________________________________
%  Do not modify anything below this line
%_____________________________________________

% At some point someone may want to change the size of the clusters used in
% the composite plots -- especially if they start wanting to do releases in
% the estuary where the shoreline is pretty tight. The following numbers
% will be the radius of the cluster in meters (origins of particle pathways
% that will be colored red), and the radius of the smaller cluster (origins
% of particle pathways colored blue).
clustsz.big = 500;
clustsz.small = 150;

switch 1
    case size(relshape,1)==1
        plottype = 'Simulation'
        reltime = reltime(1,:);
        partnum = 1;
        figure_title = ['Projected Drifter Trajectory For ' num2str(reltime(1,2)) '/' num2str(reltime(1,3)) ' Deployment'];
    case size(relshape,1)==2
        plottype = 'Transect'
        figure_title = ['Potential Drifter Trajectories From Transect Positions (' num2str(reltime(1,2)) '/' num2str(reltime(1,3)) ' Deployment)'];
        reltime = reltime(1,:);
    case size(relshape,1)>=3;
        plottype = 'Compound'
        figure_title = ['Potential Drifter Trajectories By Location & Time (' num2str(reltime(1,2)) '/' num2str(reltime(1,3)) ' Deployment)'];
end

% Determine solar-day and write to string
solday = num2str(date2jd(reltime(1,1),reltime(1,2),reltime(1,3))-date2jd(reltime(1,1)));
solday_p1 = num2str(date2jd(reltime(1,1),reltime(1,2),reltime(1,3))-date2jd(reltime(1,1))+1);

% % Create symbolic links to the input data-files 
eval(['!ln -sf /disk/ambfs19/0/workspace/forecasts/f22/' num2str(reltime(1,1)) '-' solday '/run/hgrid.gr3'])
eval(['!ln -sf /disk/ambfs19/0/workspace/forecasts/f22/' num2str(reltime(1,1)) '-' solday '/run/hgrid.ll'])
eval(['!ln -sf /disk/ambfs19/0/workspace/forecasts/f22/' num2str(reltime(1,1)) '-' solday '/run/vgrid.in'])
eval(['!ln -sf /disk/ambfs19/0/workspace/forecasts/f22/' num2str(reltime(1,1)) '-' solday '/run/param.in'])
eval(['!ln -sf /disk/ambfs19/0/workspace/forecasts/f22/' num2str(reltime(1,1)) '-' solday '/run/1_elev.61'])
eval(['!ln -sf /disk/ambfs19/0/workspace/forecasts/f22/' num2str(reltime(1,1)) '-' solday '/run/1_hvel.64'])
eval(['!ln -sf /disk/ambfs19/0/workspace/forecasts/f22/' num2str(reltime(1,1)) '-' solday '/run/1_vert.63'])
eval(['!ln -sf /disk/ambfs19/0/workspace/forecasts/f22/' num2str(reltime(1,1)) '-' solday '/run/2_elev.61'])
eval(['!ln -sf /disk/ambfs19/0/workspace/forecasts/f22/' num2str(reltime(1,1)) '-' solday '/run/2_hvel.64'])
eval(['!ln -sf /disk/ambfs19/0/workspace/forecasts/f22/' num2str(reltime(1,1)) '-' solday '/run/2_vert.63'])
eval(['!ln -sf /disk/ambfs19/0/workspace/forecasts/f22/' num2str(reltime(1,1)) '-' solday '/run/3_elev.61'])
eval(['!ln -sf /disk/ambfs19/0/workspace/forecasts/f22/' num2str(reltime(1,1)) '-' solday '/run/3_hvel.64'])
eval(['!ln -sf /disk/ambfs19/0/workspace/forecasts/f22/' num2str(reltime(1,1)) '-' solday '/run/3_vert.63'])
% disp(['!ln -sf /disk/ambfs19/0/workspace/forecasts/f22/'/' num2str(reltime(1,1)) '-' solday_p1 '/run/2_elev.61 3_elev.61']);
% eval(['!ln -sf /disk/ambfs19/0/workspace/forecasts/f22/'/' num2str(reltime(1,1)) '-' solday_p1 '/run/2_elev.61 3_elev.61'])
% eval(['!ln -sf /disk/ambfs19/0/workspace/forecasts/f22/'/' num2str(reltime(1,1)) '-' solday_p1 '/run/2_hvel.64 3_hvel.64'])
% eval(['!ln -sf /disk/ambfs19/0/workspace/forecasts/f22/'/' num2str(reltime(1,1)) '-' solday_p1 '/run/2_vert.63 3_vert.63'])

% Read the model parameterizations important to particle tracking
modpar = ptparamread('param.in');

% Process the gridfiles
fg = hgrid2fg('hgrid.gr3');
fgll = hgrid2fg('hgrid.ll');
fgll.bndpth = fg.bndpth; fgll.bnd = fg.bnd;

% Convert the date information to seconds
relsec = reltime(:,4).*3600 + reltime(:,5).*60;

% Check for lat/lons to convert to x/y coordinates
if relshape(1,1)<=90
    if size(relshape,2)==4
        relshape = [relshape(:,1)+relshape(:,2)./60 -(abs(relshape(:,3))+relshape(:,4)./60)];
    end
    for i = 1:size(relshape,1)
        dist = sqrt((relshape(i,2)-fgll.x).^2 + (relshape(i,1)-fgll.y).^2);
        [d,ind] = sort(dist);
        d = d(1:4);
        ind = ind(1:4);
        xl = min(fgll.x(ind));
        xli = find(fgll.x==xl);xli=xli(1);
        xh = max(fgll.x(ind));
        xhi = find(fgll.x==xh);xhi=xhi(1);
        xrat = (relshape(i,2)-fgll.x(xli)) / (fgll.x(xhi)-fgll.x(xli));
        yl = min(fgll.y(ind));
        yli = find(fgll.y==yl);yli = yli(1);
        yh = max(fgll.y(ind));
        yhi = find(fgll.y==yh);yhi=yhi(1);
        yrat = (relshape(i,1)-fgll.y(yli)) / (fgll.y(yhi)-fgll.y(yli));
        relshape(i,:) = [(fg.x(xli) + xrat*(fg.x(xhi)-fg.x(xli))) (fg.y(yli) + yrat*(fg.y(yhi)-fg.y(yli)))];
    end
end

% Generate the initial positions of the particles and write to a bp-file
partpos = layparticles(fg,modpar,partnum,relshape,relsec,clustsz);

% Run the ptrack_selfe application
%!./ptrack

% Wait until the run has completed
pause(420);
working = true;
while working
    fid = fopen('status.txt','r');
    out = fgets(fid);
    if ~isempty(out)
        working=false;
    end
    fclose(fid);
    pause(5);
end

% Clean house
!rm *.6* status.txt *.in

% Plot the results
hpath = plotpath(plottype,relshape,figure_title,clustsz);
sdat = hgexport('factorystyle');
hgexport(gcf,'ptrack_plot.eps',sdat,'format','eps');
!rm *.gr3 *.ll fort.*
