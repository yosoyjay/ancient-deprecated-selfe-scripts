function hpath = parttrack(modname,partnum,relshape,reltime)
% 
% h = PARTTRACK(modelname,particles,release_shape,release_time)
%
% Matlab function for generating and plotting particle tracking forecasts
% (all necessary m-files are assumed to be available in the same file as this
% function's m-file). In the case of positions, data may be supplied in either
% Lat/Lon (decimal degrees or deg decmin) or x/y pairs. This function was
% developed and tested to run on FCAST00.
%
% This is the name of the forecast as it appears in the directory call ('dev', 'db16', etc.)
% example: modelname='dev';
%
% The approximate number of particles desired per release. If relshape only
% includes one point, this value will default to 1. If relshape includes 2
% points (a transect), this number defines how many evenly spaced release
% locations there are along the transect. If relshape includes 3+
% positions, partnum defines how many particles will be released as a
% cluster from each release position
%(total particles = partnum*releasepositions*repetitions)
% example: particles=40;
%
% Defines the region in which particles are released. If only 1 point is
% listed, the particle will be released there; 2 points makes a release transect; 
% 3+ points will create a composite plot showing the projected trajectories for each
% point. NOTE: When creating composite plots, avoid choosing release points that
% lie withing 500 meters of the coastline.
% example: release_shape=[46 16 -124 5.9733
%                         46 17 -124 5.9733];
%
% When to release particles (year, month, day, hour (PST), min . Multiple releases
% are created by adding new rows with the same format. Currently is not set up to
% span 2 forecast days. Note: reltime should only include multiple release
% times when building composite plots (i.e., relshape includes 3+ positions)
% example: release_time=[2008 9 1 7 00];
%
% If a single path is requested, an ascii text file named drifterpath.dat will 
% created containing the times and positions of the particle over the course
% of the run.
% 
% In all cases, a png image of the plot will be created
%  
% Until ptparamread.m gets updated, this code only runs on the dev forecast
%

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
        plottype = 'Simulation';
        reltime = reltime(1,:);
        partnum = 1;
        figure_title = ['Projected Drifter Trajectory For ' num2str(reltime(1,2)) '/' num2str(reltime(1,3)) ' Deployment'];
    case size(relshape,1)==2
        plottype = 'Transect';
        figure_title = ['Potential Drifter Trajectories From Transect Positions (' num2str(reltime(1,2)) '/' num2str(reltime(1,3)) ' Deployment)'];
        reltime = reltime(1,:);
    case size(relshape,1)>=3;
        plottype = 'Compound';
        figure_title = ['Potential Drifter Trajectories By Location & Time (' num2str(reltime(1,2)) '/' num2str(reltime(1,3)) ' Deployment)'];
end

% Determine solar-day and write to string
solday = num2str(date2jd(reltime(1,1),reltime(1,2),reltime(1,3))-date2jd(reltime(1,1)));

% Create symbolic links to the input data-files 
eval(['!ln -s /home/workspace/local0/forecasts/' modname '/' num2str(reltime(1,1)) '-' solday '/run/hgrid.gr3'])
eval(['!ln -s /home/workspace/local0/forecasts/' modname '/' num2str(reltime(1,1)) '-' solday '/run/hgrid.ll'])
eval(['!ln -s /home/workspace/local0/forecasts/' modname '/' num2str(reltime(1,1)) '-' solday '/run/vgrid.in'])
eval(['!ln -s /home/workspace/local0/forecasts/' modname '/' num2str(reltime(1,1)) '-' solday '/run/param.in'])
eval(['!ln -s /home/workspace/local0/forecasts/' modname '/' num2str(reltime(1,1)) '-' solday '/run/1_elev.61'])
eval(['!ln -s /home/workspace/local0/forecasts/' modname '/' num2str(reltime(1,1)) '-' solday '/run/1_hvel.64'])
eval(['!ln -s /home/workspace/local0/forecasts/' modname '/' num2str(reltime(1,1)) '-' solday '/run/1_vert.63'])
eval(['!ln -s /home/workspace/local0/forecasts/' modname '/' num2str(reltime(1,1)) '-' solday '/run/2_elev.61'])
eval(['!ln -s /home/workspace/local0/forecasts/' modname '/' num2str(reltime(1,1)) '-' solday '/run/2_hvel.64'])
eval(['!ln -s /home/workspace/local0/forecasts/' modname '/' num2str(reltime(1,1)) '-' solday '/run/2_vert.63'])
eval(['!ln -s /home/workspace/local0/forecasts/' modname '/' num2str(reltime(1,1)) '-' solday '/run/3_elev.61'])
eval(['!ln -s /home/workspace/local0/forecasts/' modname '/' num2str(reltime(1,1)) '-' solday '/run/3_hvel.64'])
eval(['!ln -s /home/workspace/local0/forecasts/' modname '/' num2str(reltime(1,1)) '-' solday '/run/3_vert.63'])

% Read the model parameterizations important to particle tracking
modpar = ptparamread('param.in');

% Process the gridfiles
fg = hgrid2fg('hgrid.gr3');
fgll = hgrid2fg('hgrid.ll');
fgll.bndpth = fg.bndpth; fgll.bnd = fg.bnd;

% Convert the date information to seconds
relsec = reltime(:,4).*3600 + reltime(:,5).*60;

% Check for lat/lons to convert
switch 1
    case relshape(1,1)<=90 & size(relshape,2)==2
    % Convert decimal degrees to x/y
    for i = 1:size(relshape,1)
        dist = sqrt((relshape(i,2)-fgll.x).^2 + (relshape(i,1)-fgll.y).^2);
        [d,in] = sort(dist);
        d = d(1:3);
        ind = in(1:3);
        rat = (1./d)./sum(1./d);
        relshape(i,:) = [sum(fg.x(ind).*rat) sum(fg.y(ind).*rat)];
        aux = relshape;
    end
    case relshape(1,1)<=90 && size(relshape,2)==4
    % Convert degrees+decmin to xy
    relshape = [relshape(:,1)+relshape(:,2)./60 -(abs(relshape(:,3))+relshape(:,4)./60)];
    for i = 1:size(relshape,1)
        dist = sqrt((relshape(i,2)-fgll.x).^2 + (relshape(i,1)-fgll.y).^2);
        [d,in] = sort(dist);
        d = d(1:3);
        ind = in(1:3);
        rat = (1./d)./sum(1./d);
        relshape(i,:) = [sum(fg.x(ind).*rat) sum(fg.y(ind).*rat)];
        aux = relshape;
    end
end
if strcmp(plottype,'Compound') && typedata(1,1)<=90
    switch 1
        case relshape(1,1)<=90 & size(relshape,2)==4
        % Convert degrees+decmin to xy
        relshape = [relshape(:,1)+relshape(:,2)./60 -(abs(relshape(:,3))+relshape(:,4)./60)];
        for i = 1:size(relshape,1)
            dist = sqrt((relshape(i,2)-fgll.x).^2 + (relshape(i,1)-fgll.y).^2);
            [d,in] = sort(dist);
            d = d(1:3);
            ind = in(1:3);
            rat = (1./d)./sum(1./d);
            relshape(i,:) = [sum(fg.x(ind).*rat) sum(fg.y(ind).*rat)];
            aux = relshape;
        end
        case relshape(1,1)<=90 & size(relshape,2)==2
        % Convert decimal degrees to x/y
        for i = 1:size(relshape,1)
            dist = sqrt((relshape(i,2)-fgll.x).^2 + (relshape(i,1)-fgll.y).^2);
            [d,in] = sort(dist);
            d = d(1:3);
            ind = in(1:3);
            rat = (1./d)./sum(1./d);
            relshape(i,:) = [sum(fg.x(ind).*rat) sum(fg.y(ind).*rat)];
            aux = relshape;
        end
    end
end

% Generate the initial positions of the particles and write to a bp-file
partpos = layparticles(fg,modpar,partnum,relshape,relsec,clustsz);

% Run the ptrack_selfe application
!./ptrack

% Wait until the run has completed
pause(120);
working = true;
while working
    fid = fopen('status.txt','r');
    out = strtrim(fgets(fid));
    if strcmp(out,'Completed')
        working=false;
    end
    fclose(fid);
    pause(5);
end

% Clean house
!rm *.6* status.txt *.in

% Plot the results
hpath = plotpath(plottype,aux,figure_title);
sdat = hgexport('factorystyle');
hgexport(gcf,'ptrack_plot.eps',sdat,'format','eps');
!rm *.gr3 *.ll fort.*
