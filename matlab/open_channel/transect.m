function [] = plotTransect(hgridPath, bpPath, elevPath, varPath, timeSteps)
% This function creates a plot for a transect based on the example script in svn
%
% Input: 
%   hgridPath - Path to hgrid.gr3
%   bpPath    - Path to build point (.bp) file defining the transect
%   elevPath  - Path to ?_elev.61 file for the day of plotting     
%   varPath   - Path to the file to be plotted.  
%   timeSteps - A vector of the set of time steps to be plotted
%
% Output:
%   hPlots    - Handle to plots
%
% lopezj - 12/01/11

addpath('/usr/local/cmop/matlab/cmop/m-elio');

% Build points file
tr.hgrid = gr_readHGrid(bpPath);

% Deal with h- and v-grid
gr.hgrid = gr_readHGrid(hgridPath);
hVar  = sz_readHeader(varPath);
hElev = sz_readHeader(elevPath);
gr.vgrid = hVar.vgrid;

% Compute transect
[ob]= ob_ini_fromTrasect(gr, bpPath);
trLen = cumsum(sqrt((ob.xy.x(2:end)-ob.xy.x(1:end-1)).^2 + ...
                    (ob.xy.y(2:end)-ob.xy.y(1:end-1)).^2));
trLen = [0; trLen];

% Create a plot for every 4 time steps
for i=1:1:24
    % Read timestep of data
    [varData varTs] = sz_readTimeStep(hVar,i);
    u = varData(:,1);
    v = varData(:,2);

    % Map sz levels to depths
    du = map_sz2hts(hVar, u);
    dv = map_sz2hts(hVar, v);

    % Do something and you just end up with the data that you want in transect?
    u = ob.xy.H*double(du);
    %v = ob.xy.H*double(dv);
    v = zeros(size(u));

    % Construct vertical grid
    % Read timestep of elevations
    elevData = sz_readTimeStep(hElev,i);
    elev     = ob.xy.H*double(elevData);
    depths   = ob.xy.H*gr.hgrid.depth;
    sz=sz_computeZlevels(depths,elev,gr.vgrid);

    % x values required for every node level
    x = repmat(ob.xy.x, 1, size(sz,2)); 

    % Plotting
    figH = figure;
    plot(x,sz,'color',[0.7 0.7 0.7],'linewidth', 0.1); hold on;
    plot(x(:,end), sz(:,end), 'b', 'linewidth', 0.1)

    % Create quivers 
    plotH = quiver(x(2:5,2:2:end), sz(2:5,2:2:end),    ...
                   u(2:5,2:2:end)*100, v(2:5,2:2:end), ...
                   0,                                  ...
                   'MaxHeadSize', 0.01,                ...
                   'color', 'r');

    axis([0 1500 -10 20]);
    buf = sprintf('Time step %d', i);
    title(buf);hold on;
    % cb = colorbar('EastOutside');
    % caxis([5 20]); hold on;

    iw = 1024;
    ih = 800;

    fn = sprintf('%s-%s-%08d.png', 'oc', 'hvel', i * 900);

    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 iw/100.0 ih/100.0])
    print('-dpng', fn, '-r100');
    close(figH);
end
