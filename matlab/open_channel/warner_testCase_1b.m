function [] = warner_testCase_1b(runPath, movie)
% Creates plots similar to those found in Warner et al. 2005 for test case 1b. 
%
%
% Input:     
%   runPath - Path to run directory, assumes outputs are already combined in outputs dir
%           - Assumes *.bp file named oc_length.bp in directory
%   movie   - 1 - Makes a avi, 0 - No avi
%
% Output:
%
% lopezj - 12/08/11

% Constants etc.
tic
addpath('/usr/local/cmop/matlab/cmop/m-elio');
xPos = 4;           % Index of x position in open channel to get data
H = 10;             % Total depth of water column
botRough = 0.0005   % Bottom roughness 
gamma = 2e-4;       % Slope
anUBar = 1;         % Depth averaged horizontal velocity 
%Cd = 1.e-2;        % Drag coefficient
%grav = 9.81        % Gravity 
% dahv = sqrt(grav*gamma*H/Cd);   % Depth average velocity
% lambda = grav*gamma/dahv*(log(H/z1)-1);

% Build points file
bpPath = sprintf('%s/oc_length.bp', runPath);
tr.hgrid = gr_readHGrid(bpPath);

% Deal with h- and v-grid
% Note that sz layers must be extracted from binary output file headers
% Open up files for all variables I'm looking at in open channel
hgridPath = sprintf('%s/hgrid.gr3', runPath);
gr.hgrid = gr_readHGrid(hgridPath);
hvelPath = sprintf('%s/outputs/1_hvel.64', runPath);
hHvel = sz_readHeader(hvelPath);
kinePath = sprintf('%s/outputs/1_kine.63', runPath);
hKine = sz_readHeader(kinePath);
tdffPath = sprintf('%s/outputs/1_tdff.63', runPath);
hTdff = sz_readHeader(tdffPath);
vdffPath = sprintf('%s/outputs/1_vdff.63', runPath);
hVdff = sz_readHeader(vdffPath);
elevPath = sprintf('%s/outputs/1_elev.61', runPath);
hElev = sz_readHeader(elevPath);
gr.vgrid = hKine.vgrid;

% Compute transect
[ob]= ob_ini_fromTrasect(gr, bpPath);
trLen = cumsum(sqrt((ob.xy.x(2:end)-ob.xy.x(1:end-1)).^2 + ...
                    (ob.xy.y(2:end)-ob.xy.y(1:end-1)).^2));
trLen = [0; trLen];

% Create animation 
if movie == 1
    movieH = avifile('open_channel.avi');
end

% Create a plot for every 4 time steps
for i=4:4:96
    % Read timestep of data
    [hvelData varTs] = sz_readTimeStep(hHvel,i);
    [kineData varTs] = sz_readTimeStep(hKine,i);
    [tdffData varTs] = sz_readTimeStep(hTdff,i);
    [vdffData varTs] = sz_readTimeStep(hVdff,i);

    % Map sz levels to depths
    dHvel = map_sz2hts(hHvel, hvelData(:,1));
    dKine = map_sz2hts(hKine, kineData);
    dTdff = map_sz2hts(hTdff, tdffData);
    dVdff = map_sz2hts(hVdff, vdffData);

    % Do this and you just end up with the data that you want in transect
    % ob is the object that has the transect information
    hvel = ob.xy.H*double(dHvel);
    kine = ob.xy.H*double(dKine);
    tdff = ob.xy.H*double(dTdff);
    vdff = ob.xy.H*double(dVdff);

    % Construct vertical grid
    % Read timestep of elevations
    elevData = sz_readTimeStep(hElev,i);
    elev     = ob.xy.H*double(elevData);
    depths   = ob.xy.H*gr.hgrid.depth;
    sz=sz_computeZlevels(depths,elev,gr.vgrid);

    % Copy x values so every node vertical level has an x value
    x = repmat(ob.xy.x, 1, size(sz,2)); 


    % Calculate analytical solution via Warner et al. [2005]
    % All depths sz are converted to 0 at surface reference to at bottom
    % tdff - Eddy diffusivity
    anUStar = (0.41*anUBar) / ...               
              (log(H/botRough)-1+(botRough/H));
    a = 0.41*anUStar*(H+sz);
    b = (1-((H+sz)./H));
    anTdff = (a.*b)/0.8;
    % vdff - Eddy viscosity
    a = 0.41*0.061*(H+sz);
    b = (1-((H+sz)./H));
    anVdff = (a.*b);
    % hvel - Horizontal velocity - Only U component as the other is so small 
    % All analaytical calcs from Warner et al. [2005] open channel test
    anUStar = (0.41*anUBar) / ...               
              (log(H/botRough)-1+(botRough/H));
    % anSol is the analytical velocity
    anHvel = ((1/0.41)*log((H+sz)/botRough))*anUStar;
    anHvel(isinf(anHvel)) = 0;
    % Bottom level ends up in inf, so I replace it with zero which it replaces
    % kine - Turbulent kinetic energy

    % Axis limits for all variables
    axisDimsTdff = [0 0.08 -10 1];
    axisDimsVdff = [0 0.08 -10 1];
    axisDimsHvel = [0 1.5 -10 1];
    axisDimsKine = [0 0.015 -10 1];

    % Plot data and analytical solution
    % 4 sub plots for each variable
    screen_sz = get(0,'ScreenSize');
    figH = figure('Position',[ 1 screen_sz(4)/2 1200 800], 'visible', 'off');
   
    % Velocity
    h(1) = subplot(1,4,1);
    sp_pos = get(h(1), 'position'); hold on;
    ph(1)   = plot(-hvel(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color','k', 'marker','*'); 
    anPh(1) = plot(anHvel(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color', 'r', ...
                   'marker', 'o', 'linestyle', '--');
    % Calculate RMSE and add to plot
    % u is neg vel. anSol positive so the + below is correct
    err = hvel(xPos,:)+anHvel(xPos,:);
    mse = sum(err.^2/numel(hvel(xPos,:)));
    rmse = sqrt(mse);
    % Clean up 
    axis(axisDimsHvel);
    buf = sprintf('RMSE %f', rmse);
    title(buf);
    %legH = legend('Model', 'Analytical', 'Location', 'NorthEast');
    ylabel('Depth m');
    xlabel('Velocity m/s');

    % Plot levels for reference
    plot(x,sz,'color',[0.7 0.7 0.7],'linewidth', 0.1); hold on;
    plot(x(:,end), sz(:,end), 'b', 'linewidth', 0.1);

    % TKE
    h(2) = subplot(1,4,2);
    sp_pos = get(h(2), 'position'); hold on;
    ph(2)   = plot(kine(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color','k', 'marker','*'); 
    set(gca,'XMinorTick', 'on');
    %anPh(2) = plot(anKine(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color', 'r', ...
    %               'marker', 'o', 'linestyle', '--');
    % Calculate RMSE and add to plot
    %err = kine(xPos,:)-anKine(xPos,:);
    %mse = sum(err.^2/numel(kine(xPos,:)));
    %rmse = sqrt(mse);
    % Clean up 
    axis(axisDimsKine);
    %buf = sprintf('TKE - RMSE %f', i, rmse);
    %title(buf);
    %legH = legend('Model', 'Analytical', 'Location', 'NorthEast');
    ylabel('Depth m');
    xlabel('TKE m^2/s^2');
    % Plot levels for reference
    plot(x,sz,'color',[0.7 0.7 0.7],'linewidth', 0.1); hold on;
    plot(x(:,end), sz(:,end), 'b', 'linewidth', 0.1);

    % Eddy diffusivity 
    h(3) = subplot(1,4,3);
    sp_pos = get(h(3), 'position'); hold on;
    ph(3)   = plot(tdff(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color','k', 'marker','*'); 
    anPh(3) = plot(anTdff(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color', 'r', ...
                   'marker', 'o', 'linestyle', '--');
    % Calculate RMSE and add to plot
    err = tdff(xPos,:)-anTdff(xPos,:);
    mse = sum(err.^2/numel(tdff(xPos,:)));
    rmse = sqrt(mse);
    % Clean up 
    axis(axisDimsTdff);
    buf = sprintf('RMSE %f', rmse);
    title(buf);
    %legH = legend('Model', 'Analytical', 'Location', 'NorthEast');
    ylabel('Depth m');
    xlabel('Eddy diffusivity, m^2/s');
    % Plot levels for reference
    plot(x,sz,'color',[0.7 0.7 0.7],'linewidth', 0.1); hold on;
    plot(x(:,end), sz(:,end), 'b', 'linewidth', 0.1);

    % Eddy viscosity 
    h(4) = subplot(1,4,4);
    sp_pos = get(h(4), 'position'); hold on;
    ph(4)   = plot(vdff(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color','k', 'marker','*'); 
    anPh(4) = plot(anVdff(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color', 'r', ...
                   'marker', 'o', 'linestyle', '--');
    % Calculate RMSE and add to plot
    err = vdff(xPos,:)-anVdff(xPos,:);
    mse = sum(err.^2/numel(vdff(xPos,:)));
    rmse = sqrt(mse);
    % Clean up 
    axis(axisDimsVdff);
    buf = sprintf('RMSE %f', rmse);
    title(buf);
    %legH = legend('Model', 'Analytical', 'Location', 'NorthEast');
    ylabel('Depth m');
    xlabel('Eddy viscosity, m^2/s');
    % Plot levels for reference
    plot(x,sz,'color',[0.7 0.7 0.7],'linewidth', 0.1); hold on;
    plot(x(:,end), sz(:,end), 'b', 'linewidth', 0.1);

    % Get ready for saves
    fn = sprintf('%s-%s-%06d.png', 'oc', 'all_vars', i * 900);
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
toc
