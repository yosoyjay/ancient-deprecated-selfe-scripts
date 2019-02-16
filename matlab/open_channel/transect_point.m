function [] = plotTransect(hgridPath, bpPath, elevPath, varPath, timeSteps, movie)
% This function creates a plot for a transect based on the example script in svn
%
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

% Constants etc.
addpath('/usr/local/cmop/matlab/cmop/m-elio');
xPos = 4;       	% Index of x position in open channel to get data
H = 10;         	% Total depth of water column
botRough = 0.00053 	% Bottom roughness 
gamma = 2e-4;   	% Slope
anUBar = 1;        % Depth averaged horizontal velocity 
Cd = 1.e-2;     	% Drag coefficient
grav = 9.81     	% Gravity 
% dahv = sqrt(grav*gamma*H/Cd);   % Depth average velocity
% lambda = grav*gamma/dahv*(log(H/z1)-1);

% Build points file
tr.hgrid = gr_readHGrid(bpPath);

% Deal with h- and v-grid
% Note that sz layers must be extracted from binary output file headers
gr.hgrid = gr_readHGrid(hgridPath);
hVar  = sz_readHeader(varPath);
hElev = sz_readHeader(elevPath);
gr.vgrid = hVar.vgrid;

% Compute transect
[ob]= ob_ini_fromTrasect(gr, bpPath);
trLen = cumsum(sqrt((ob.xy.x(2:end)-ob.xy.x(1:end-1)).^2 + ...
                    (ob.xy.y(2:end)-ob.xy.y(1:end-1)).^2));
trLen = [0; trLen];

% Determine plot variable and appropriate lables, axis, etc...
varSPos = strfind(varPath, '_');
varName = varPath(varSPos+1:varSPos+4);
if varName == 'salt'
    dim = 1;
    axisDims = [0 32 -10 1];
elseif varName == 'temp'
    dim = 1;
    axisDims = [5 10 -10 1];
elseif varName == 'tdff'
    dim = 1;
    axisDims = [0 0.08 -10 1];
elseif varName == 'vdff'
    dim = 1;
    axisDims = [0 0.08 -10 1];
elseif varName == 'hvel'
    dim = 2;
    axisDims = [0 1.5 -10 1];
elseif varName == 'dahv'
    dim = 2;
    axisDims = [-3 3 -10 1];
elseif varName == 'kine'
    dim = 1;
    axisDims = [0 0.1 -10 1];
else
    dim = 1;
    axisDims = [0 0.25 -10 1];
end

% Create animation 
if movie == 1
    movieF = sprintf('%s.avi', varName); 
    movieH = avifile(movieF);
end

% Create a plot for every 4 time steps
for i=4:4:40
    % Read timestep of data
    [varData varTs] = sz_readTimeStep(hVar,i);
    if dim == 2
        u = varData(:,1);
        v = varData(:,2);
    else
        u = varData(:,1);
    end

    % Map sz levels to depths
    du = map_sz2hts(hVar, u);

    % Do something and you just end up with the data that you want in transect?
    u = ob.xy.H*double(du);
    v = zeros(size(u));

    % Construct vertical grid
    % Read timestep of elevations
    elevData = sz_readTimeStep(hElev,i);
    elev     = ob.xy.H*double(elevData);
    depths   = ob.xy.H*gr.hgrid.depth;
    sz=sz_computeZlevels(depths,elev,gr.vgrid);

    % x values required for every node level
    x = repmat(ob.xy.x, 1, size(sz,2)); 

    % Just get values for 800 m from boundary 
    % u = u(3,:);
    % sz = sz(3,:);

    % Calculate analytical solution via Warner et al. [2005]
    % All depths sz are converted to 0 at surface reference to at bottom
    if varName == 'salt'
    elseif varName == 'temp'
    elseif varName == 'tdff'
        anUStar = (0.41*anUBar) / ...               
                  (log(H/botRough)-1+(botRough/H));
        a = 0.41*anUStar*(H+sz);
        b = (1-((H+sz)./H));
        anSol = (a.*b)/0.8;
    elseif varName == 'vdff'
        a = 0.41*0.061*(H+sz);
        b = (1-((H+sz)./H));
        anSol = (a.*b);
    elseif varName == 'hvel'
        % All analaytical calcs from Warner et al. [2005] open channel test
        anUStar = (0.41*anUBar) / ...               
                  (log(H/botRough)-1+(botRough/H));
        % anSol is the analytical velocity
        anSol = ((1/0.41)*log((H+sz)/botRough))*anUStar;
        % Bottom level ends up in inf, so I replace it with zero which it replaces
        anSol(isinf(anSol)) = 0;
    elseif varName == 'dahv'
        dim = 1;
        axisDims = [0 3 -10 1];
        anSol = nan(size(sz))*-10;
    elseif varName == 'kine'
        dim = 1;
        axisDims = [0 0.1 -10 1];
        anSol = nan(size(sz))*-10;
    else
        dim = 1;
        axisDims = [0 0.25 -10 1];
    end

    % Plot data and analytical solution
    if dim == 2
        plotH = plot(-u(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color','k', 'marker','*'); 
        panaH = plot(anSol(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color', 'r', ...
                     'marker', 'o', 'linestyle', '--');
        % Calculate RMSE and add to plot
        % u is neg vel. anSol positive so the + below is correct
        err = u(xPos,:)+anSol(xPos,:);
        mse = sum(err.^2/numel(u(xPos,:)));
        rmse = sqrt(mse);

    else
		u(xPos,:)
        plotH = plot(u(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color','k', 'marker','*'); 
        panaH = plot(anSol(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color', 'r', ...
                     'marker', 'o', 'linestyle', '--');
        % Calculate RMSE and add to plot
        err = u(xPos,:)-anSol(xPos,:);
        mse = sum(err.^2/numel(u(xPos,:)));
        rmse = sqrt(mse);

    end
   
    axis(axisDims);
    buf = sprintf('Time step %d: RMSE %f', i, rmse);
    title(buf);hold on;
    legH = legend('Model', 'Analytical', 'Location', 'NorthWest');

    % Plot levels for reference
    figH = figure;
    plot(x,sz,'color',[0.7 0.7 0.7],'linewidth', 0.1); hold on;
    plot(x(:,end), sz(:,end), 'b', 'linewidth', 0.1);

    % Get ready for saves
    fn = sprintf('%s-%s-%06d.png', 'oc', varName, i * 900);
    iw = 1024;
    ih = 800;
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 iw/100.0 ih/100.0])
    print('-dpng', fn, '-r100');
    if movie == 1
        movieH = addframe(movieH, figH);
    end
    close(figH);
end
if movie == 1
    movieH = close(movieH);
end
