function [] = plotTransect(runPath, movie)
% This function creates a plot for a transect based on the example script in svn
%
%
% Input:	 
%	runPath - Path to run directory, assumes outputs are already combined in outputs dir
%			- Assumes *.bp file named oc_length.bp in directory
%	movie 	- 1 - Makes a avi, 0 - No avi
%
% Output:
%
% lopezj - 12/08/11

% Constants etc.
addpath('/usr/local/cmop/matlab/cmop/m-elio');
xPos = 4;       	% Index of x position in open channel to get data
H = 10;         	% Total depth of water column
botRough = 0.00053 	% Bottom roughness 
gamma = 2e-4;   	% Slope
anUBar = 1;         % Depth averaged horizontal velocity 
%Cd = 1.e-2;     	% Drag coefficient
%grav = 9.81     	% Gravity 
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

% KL
hvelPathKL = sprintf('%s/outputs_kl/1_hvel.64', runPath);
hHvelKL = sz_readHeader(hvelPathKL);
kinePathKL = sprintf('%s/outputs_kl/1_kine.63', runPath);
hKineKL = sz_readHeader(kinePathKL);
tdffPathKL = sprintf('%s/outputs_kl/1_tdff.63', runPath);
hTdffKL = sz_readHeader(tdffPathKL);
vdffPathKL = sprintf('%s/outputs_kl/1_vdff.63', runPath);
hVdffKL = sz_readHeader(vdffPathKL);
elevPathKL = sprintf('%s/outputs_kl/1_elev.61', runPath);
hElevKL = sz_readHeader(elevPathKL);
sedPathKL = sprintf('%s/outputs_kl/1_trcr_1.63', runPath);
hSedKL = sz_readHeader(sedPathKL);

% E 
hvelPathE = sprintf('%s/outputs_e/1_hvel.64', runPath);
hHvelE = sz_readHeader(hvelPathE);
kinePathE = sprintf('%s/outputs_e/1_kine.63', runPath);
hKineE = sz_readHeader(kinePathE);
tdffPathE = sprintf('%s/outputs_e/1_tdff.63', runPath);
hTdffE = sz_readHeader(tdffPathE);
vdffPathE = sprintf('%s/outputs_e/1_vdff.63', runPath);
hVdffE = sz_readHeader(vdffPathE);
elevPathE = sprintf('%s/outputs_e/1_elev.61', runPath);
hElevE = sz_readHeader(elevPathE);
sedPathE = sprintf('%s/outputs_e/1_trcr_1.63', runPath);
hSedE = sz_readHeader(sedPathE);

% W 
hvelPathW = sprintf('%s/outputs_w/1_hvel.64', runPath);
hHvelW = sz_readHeader(hvelPathW);
kinePathW = sprintf('%s/outputs_w/1_kine.63', runPath);
hKineW = sz_readHeader(kinePathW);
tdffPathW = sprintf('%s/outputs_w/1_tdff.63', runPath);
hTdffW = sz_readHeader(tdffPathW);
vdffPathW = sprintf('%s/outputs_w/1_vdff.63', runPath);
hVdffW = sz_readHeader(vdffPathW);
elevPathW = sprintf('%s/outputs_w/1_elev.61', runPath);
hWlevE = sz_readHeader(elevPathW);
sedPathW = sprintf('%s/outputs_w/1_trcr_1.63', runPath);
hSedW = sz_readHeader(sedPathE);

% Grab vgrid from one file
gr.vgrid = hKineKL.vgrid;


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
for i=4:4:40
    % Read timestep of data
	% KL
    [hvelKLData varTs] = sz_readTimeStep(hHvelKL,i);
    [kineKLData varTs] = sz_readTimeStep(hKineKL,i);
    [tdffKLData varTs] = sz_readTimeStep(hTdffKL,i);
    [vdffKLData varTs] = sz_readTimeStep(hVdffKL,i);
    [sedKLData varTs]  = sz_readTimeStep(hSedKL,i);
	% E 
    [hvelEData varTs] = sz_readTimeStep(hHvelE,i);
    [kineEData varTs] = sz_readTimeStep(hKineE,i);
    [tdffEData varTs] = sz_readTimeStep(hTdffE,i);
    [vdffEData varTs] = sz_readTimeStep(hVdffE,i);
    [sedEData varTs]  = sz_readTimeStep(hSedE,i);
	% W 
    [hvelWData varTs] = sz_readTimeStep(hHvelW,i);
    [kineWData varTs] = sz_readTimeStep(hKineW,i);
    [tdffWData varTs] = sz_readTimeStep(hTdffW,i);
    [vdffWData varTs] = sz_readTimeStep(hVdffW,i);
    [sedWData varTs]  = sz_readTimeStep(hSedW,i);

    % Map sz levels to depths
	% KL
    dHvelKL = map_sz2hts(hHvelKL, hvelKLData(:,1));
    dKineKL = map_sz2hts(hKineKL, kineKLData);
    dTdffKL = map_sz2hts(hTdffKL, tdffKLData);
    dVdffKL = map_sz2hts(hVdffKL, vdffKLData);
    dSedKL = map_sz2hts(hSedKL, sedKLData);
	% E
    dHvelE = map_sz2hts(hHvelE, hvelEData(:,1));
    dKineE = map_sz2hts(hKineE, kineEData);
    dTdffE = map_sz2hts(hTdffE, tdffEData);
    dVdffE = map_sz2hts(hVdffE, vdffEData);
    dSedE = map_sz2hts(hSedE, sedEData);
	% W
    dHvelW = map_sz2hts(hHvelW, hvelWData(:,1));
    dKineW = map_sz2hts(hKineW, kineWData);
    dTdffW = map_sz2hts(hTdffW, tdffWData);
    dVdffW = map_sz2hts(hVdffW, vdffWData);
    dSedW = map_sz2hts(hSedW, sedWData);

    % Do this and you just end up with the data that you want in transect
	% ob is the object that has the transect information
	% KL
    hvelKL = ob.xy.H*double(dHvelKL);
    kineKL = ob.xy.H*double(dKineKL);
    tdffKL = ob.xy.H*double(dTdffKL);
    vdffKL = ob.xy.H*double(dVdffKL);
    sedKL = ob.xy.H*double(dSedKL)
	% E
    hvelE = ob.xy.H*double(dHvelE);
    kineE = ob.xy.H*double(dKineE);
    tdffE = ob.xy.H*double(dTdffE);
    vdffE = ob.xy.H*double(dVdffE);
    sedE = ob.xy.H*double(dSedE)
	% W
    hvelW = ob.xy.H*double(dHvelW);
    kineW = ob.xy.H*double(dKineW);
    tdffW = ob.xy.H*double(dTdffW);
    vdffW = ob.xy.H*double(dVdffW);
    sedW = ob.xy.H*double(dSedW)

    % Construct vertical grid
    % Read timestep of elevations
    elevData = sz_readTimeStep(hElevKL,i);
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
    axisDimsKine = [0 0.05 -10 1];
	axisDimsSed  = [0 2 -10 1];

    % Plot data and analytical solution
	% 4 sub plots for each variable
    screen_sz = get(0,'ScreenSize');
 	figH = figure('Position',[ 1 screen_sz(4)/2 1200 800]);
   
	% Velocity
	h(1) = subplot(1,5,1);
	sp_pos = get(h(1), 'position'); hold on;
	ph(1)   = plot(-hvelKL(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color','r', 'marker','*'); 
	ph(2)   = plot(-hvelE(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color','b', 'marker','+'); 
	ph(3)   = plot(-hvelW(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color','g', 'marker','x'); 
    anPh(1) = plot(anHvel(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color','k', ...
                   'marker', 'o', 'linestyle', '--');
    % Calculate RMSE and add to plot
    % u is neg vel. anSol positive so the + below is correct
    %err = hvel(xPos,:)+anHvel(xPos,:);
    %mse = sum(err.^2/numel(hvel(xPos,:)));
    %rmse = sqrt(mse);
	% Clean up 
    axis(axisDimsHvel);
    %buf = sprintf('RMSE %f', rmse);
    %title(buf);
    legH = legend('k-kl', 'k-e', 'k-w', 'Analytical');
	set(legH, 'Position', [0.02 0.85 0.1 0.1]);
    ylabel('Depth m');
	xlabel('Velocity m/s');

    % Plot levels for reference
    plot(x,sz,'color',[0.7 0.7 0.7],'linewidth', 0.1); hold on;
    plot(x(:,end), sz(:,end), 'b', 'linewidth', 0.1);

	% TKE
	h(2) = subplot(1,5,2);
	sp_pos = get(h(2), 'position'); hold on;
	ph(4)   = plot(kineKL(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color','r', 'marker','*'); 
	kh(5)   = plot(kineE(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color','b', 'marker','+'); 
	ph(6)   = plot(kineW(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color','g', 'marker','x'); 
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
	xlabel('TKE m^2/s^2');
    % Plot levels for reference
    plot(x,sz,'color',[0.7 0.7 0.7],'linewidth', 0.1); hold on;
    plot(x(:,end), sz(:,end), 'b', 'linewidth', 0.1);

	% Eddy diffusivity 
	h(3) = subplot(1,5,3);
	sp_pos = get(h(3), 'position'); hold on;
	ph(7)   = plot(tdffKL(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color','r', 'marker','*'); 
	ph(8)   = plot(tdffE(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color','b', 'marker','+'); 
	ph(9)   = plot(tdffW(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color','g', 'marker','x'); 
    anPh(3) = plot(anTdff(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color','k', ...
                   'marker', 'o', 'linestyle', '--');
    % Calculate RMSE and add to plot
    %err = tdff(xPos,:)-anTdff(xPos,:);
    %mse = sum(err.^2/numel(tdff(xPos,:)));
    %rmse = sqrt(mse);
	% Clean up 
    axis(axisDimsTdff);
    %buf = sprintf('RMSE %f', rmse);
    %title(buf);
    %legH = legend('Model', 'Analytical', 'Location', 'NorthEast');
	xlabel('Eddy diffusivity, m^2/s');
    % Plot levels for reference
    plot(x,sz,'color',[0.7 0.7 0.7],'linewidth', 0.1); hold on;
    plot(x(:,end), sz(:,end), 'b', 'linewidth', 0.1);

	% Eddy viscosity 
	h(4) = subplot(1,5,4);
	sp_pos = get(h(4), 'position'); hold on;
	ph(10)   = plot(vdffKL(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color','r', 'marker','*'); 
	ph(11)   = plot(vdffE(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color','b', 'marker','+'); 
	ph(11)   = plot(vdffW(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color','g', 'marker','x'); 
    anPh(4) = plot(anVdff(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color','k', ...
                   'marker', 'o', 'linestyle', '--');
    % Calculate RMSE and add to plot
    %err = vdff(xPos,:)-anVdff(xPos,:);
    %mse = sum(err.^2/numel(vdff(xPos,:)));
    %rmse = sqrt(mse);
	% Clean up 
    axis(axisDimsVdff);
    %buf = sprintf('RMSE %f', rmse);
    %title(buf);
    %legH = legend('Model', 'Analytical', 'Location', 'NorthEast');
	xlabel('Eddy viscosity, m^2/s');
    % Plot levels for reference
    plot(x,sz,'color',[0.7 0.7 0.7],'linewidth', 0.1); hold on;
    plot(x(:,end), sz(:,end), 'b', 'linewidth', 0.1);

	% Sediment 
	h(5) = subplot(1,5,5);
	sp_pos = get(h(5), 'position'); hold on;
	ph(12)   = plot(sedKL(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color','r', 'marker','*'); 
	ph(13)   = plot(sedE(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color','b', 'marker','+'); 
	ph(14)   = plot(sedW(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color','g', 'marker','x'); 
%    anPh(5) = plot(anVdff(xPos,:), sz(xPos,:), 'linewidth', 1.5, 'color', 'r', ...
%                   'marker', 'o', 'linestyle', '--');
    % Calculate RMSE and add to plot
%    err = vdff(xPos,:)-anVdff(xPos,:);
%    mse = sum(err.^2/numel(vdff(xPos,:)));
%    rmse = sqrt(mse);
	% Clean up 
    axis(axisDimsSed);
    %buf = sprintf('RMSE %f', rmse);
%    title(buf);
    %legH = legend('Model', 'Analytical', 'Location', 'NorthEast');
	xlabel('Sediment, kg/m^3');
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
