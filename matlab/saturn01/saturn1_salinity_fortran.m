function [] = saturn1_salinity(path, start_date, n_days, model_variable, grid_num, varargin)
% [] = saturn1_salinity(start_date, n_days, model_variable, grid_num, varargin)
%
% This function makes a vertical profile from model results at Saturn1 or the N. Channel.
% Updated form previous version to use Joseph's fortran code.
%
% Input:
% path       - Path to run directory
% start_date - First day to extract data 
% n_days     - Number of days to extract
% model_variable - What variable from the model do you want? Use file names salt.63 
% grid_num   - What grid is used? (Only for plot extraction)
% varargs
%   output_dt - Number of seconds per output time step
%   cb_auto   - Turn on autoscale for colorbar  '1', else it's for salinity [0 32]  
%	official - 1 yes, 0 no
%

% Deal with varargs
opt_args = size(varargin,2);
if opt_args == 1
    output_dt = varargin{1};
    cb_auto = 0;
	just_saturn = 1;
elseif opt_args == 2
    output_dt = varargin{1};
    cb_auto = varargin{2};
	just_saturn = 1;
elseif opt_args == 3
    output_dt = varargin{1};
    cb_auto = varargin{2};
	official = varargin{3};
else
    output_dt = 960;
    cb_auto = 0;
	official = 0;
end

% Map file names to plot labels 
file_key = {'vert.63', 'temp.63', 'salt.63', 'conc.63',   ...
            'tdff.63', 'vdff.63', 'kine.63', 'mixl.63',   ...
            'zcor.63', 'trcr_1.63'};
label_val = {'Vertical Velocity (m/s)', 'Temperature (C)', ...
             'Salinity (psu)', 'Density (kg/m^3)',         ...
             'Eddy Diffusivity', 'Eddy Viscosity',         ...
             'Turbulent Kinetic Energy', 'Mixing Length'   ...
             'Z Coordinates', 'Sediment (kg/m^3)'};
labels = containers.Map(file_key, label_val);
cb_label = labels(model_variable);

% Saturn01 or A few N. Channel and mouth locations that I'm intersted in 
%if just_saturn == 0
%    v_stations = [ 333049.32, 293774.45; ...
%                   346139.87, 290622.48; ...
%                   348815.60, 290622.48; ...
%                   349556.20, 290839.20; ...
%                   350133.66, 291732.58; ...
%                   351004.42, 292215.04];
%    names = {'Mouth', 'N. Channel 1', 'N. Channel 2', 'Saturn01', ...
%             'N. Channel 3', 'N. Channel 4'};
%else 
%    names = {'Saturn01'};
%    v_stations = [ 349300, 290931 ];
%end


% Check path and availability of modext
%ret_code = exist('modext.m','file');
%if ret_code ~= 2
%    addpath '/usr/local/cmop/modext';
%    addpath '/home/workspace/users/lopezj/scripts/matlab/';
%end
%ret_code = exist('imagescnan.m','file');
%if ret_code ~= 2
%    addpath '/home/workspace/users/lopezj/scripts/matlab/plotting';
%end

path = sprintf('%s/outputs',path);
cd(path);
start_date = datenum(start_date);

% Make station.bp
fid = fopen('station.bp','w');
fprintf(fid,'Profiles in North Channel\n');
fprintf(fid,'1\n');
fprintf(fid,'1 349300 290932\n');
fclose(fid);

% Make extract.in
fid = fopen('extract.in','w');
fprintf(fid,'%s\n', model_variable);
fprintf(fid,'%d\n', n_days);
fclose(fid);

% Link to binary to extract model data
[status, result] = unix('ln -fs /home/workspace/users/lopezj/scripts/Post-Processing-Fortran/read_output7b_group ./');

%for i = 1:size(names,2)
for i = 1:1
	% Extract data and load data
	[status, result] = unix('./read_output7b_group');
	load fort.18

	% Reshape data for plot: x-time, y-depth, z-data
    nLevels = size(unique(fort(:,3)),1);
	nDT = size(fort,1)/nLevels; 

	depths = reshape(fort(:,4),nLevels,nDT);
	data = reshape(fort(:,5),nLevels,nDT);
	time = reshape(fort(:,1),nLevels,nDT);
	dates = start_date + time;

	% Plot as surface and rotate around
	fH = figure('color','w');
	surf(time,depths,data,'edgecolor','none'); 
    surf(time,depths,data,'edgecolor','none');
    hold on;
	set(gca,'layer','top');
 	grid on;
	set(gca,'fontsize',14);
	view(0,90);

	% X-Axis
	labels = [];
	for day = 0:n_days*2
        if mod(day,2) == 0
            labels = [labels; datestr(start_date + day/2, 'mm/dd')];
        else
            labels = [labels; '     '];
        end
	end
	set(gca,'xticklabel',labels);
	%datetick('x','mm/dd');

    % Y-Axis
    ylabel('Depth From MSL (m)', 'FontSize', 14);

    % Colorbar
    hcb = colorbar;
    set(get(hcb,'ylabel'),'string',cb_label, 'FontSize', 14);
    if cb_auto == 0 
        caxis([0 32]);
    end

    % Title Etc..
    grid on;
    %title_date = sprintf('%s - %s', datestr(start_date, 'mm/dd/yyyy'), ...
    %             datestr(datenum(start_date+n_days-1), 'mm/dd/yyyy'));
    % title_grid = sprintf('%s %d', names{i}, grid_num);
    title_grid = sprintf('%s %d', 'Saturn01', grid_num);
    %title({title_grid;title_date}, 'FontSize', 14);
    title({title_grid}, 'FontSize', 14);

    % Anti-alias
    % aa = myaa(4);

    % Save image as *.png to post
    %image_str = sprintf('%s_%s_%d_%s_%i.png', names{i}, model_variable, ...
    %            grid_num, datestr(start_date, 'mm-dd'), n_days);
    image_str = sprintf('%s_%s_%d_%s_%i.png', 'Saturn 01', model_variable, ...
                grid_num, datestr(start_date, 'mm-dd'), n_days);
    saveas(fH,image_str);
end
