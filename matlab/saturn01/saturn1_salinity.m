function [] = saturn1_salinity(start_date, n_days, model_variable, grid_num, varargin)
% [] = saturn1_salinity(start_date, n_days, model_variable, grid_num, varargin)
% This function makes a vertical profile from model results at Saturn1 or the N. Channel
%
% Input:
% start_date - First day to extract data 
% n_days     - Number of days to extract
% model_variable - What variable from the model do you want? Use file names salt.63 
% grid_num   - What grid is used? (Only for plot extraction)
% varargs
%   cb_auto   - Turn on autoscale for colorbar  '1', else it's for salinity [0 32]  
%   just_saturn - Plot only Saturn01 (T or F)  '1' for profile of just Saturn01 or 
%                 '0' to make a unique vertical profile for each virtual station labeled
%                 below.
%	official  - 

% Deal with varargs
opt_args = size(varargin,2);
if opt_args == 1
    cb_auto = varargin{1};
	just_saturn = 1;
	official = 0;
elseif opt_args == 2
    cb_auto = varargin{1};
    just_saturn = varargin{2};
	official = 0;
elseif opt_args == 3
    cb_auto = varargin{2};
    just_saturn = varargin{3};
	official = varargin{4};
else
    cb_auto = 0;
    just_saturn = 1;
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
if just_saturn == 0
    v_stations = [ 333049.32, 293774.45; ...
                   346139.87, 290622.48; ...
                   348815.60, 290622.48; ...
                   349556.20, 290839.20; ...
                   350133.66, 291732.58; ...
                   351004.42, 292215.04];
    names = {'Mouth', 'N. Channel 1', 'N. Channel 2', 'Saturn01', ...
             'N. Channel 3', 'N. Channel 4'};
else 
    names = {'Saturn01'};
	v_stations = [ 349300, 290931 ];
end


% Check path and availability of modext
ret_code = exist('modext.m','file');
if ret_code ~= 2
%    addpath '/usr/local/cmop/modext';
    addpath '/home/workspace/users/lopezj/scripts/matlab/';
end
ret_code = exist('imagescnan.m','file');
if ret_code ~= 2
    addpath '/home/workspace/users/lopezj/scripts/matlab/plotting';
end

% Data collected every 15 minutes at every 0.1 m to 25 meters
% n_days*output_dt - 1 because that is the actual number of available times.
start_date = datenum(start_date);
times = start_date:1/(86400/output_dt):addtodate(start_date,n_days,'day');
depths = 0:0.25:25;


for i = 1:size(names,2)
	% Grab data
	if official == 0
		salinity = modext(times(2:end-1), v_stations(i,2), v_stations(i,1), depths, ...
				   model_variable, './', output_dt);
	else
		salinity = modext(times(2:end-1), v_stations(i,2), v_stations(i,1), depths, ...
				   model_variable, '22', output_dt);
	end

	% Plot it up
	fH = figure('color','w');
    imagescnan(salinity);

    % X-Axis
    set(gca,'xtick', 0:output_dt:(n_days*output_dt));
    day_one = datestr(datenum(start_date+1),'mm/dd');
    day_two = datestr(datenum(start_date+2), 'mm/dd');
    xtick_labels = {day_two, day_one};
    set(gca,'xticklabel', xtick_labels);

    % Y-Axis
    set(gca,'ytick', 0:20:101);
    set(gca,'yticklabel', depths(1:20:101));
    ylabel('Depth From Free Surface (m)', 'FontSize', 12);

    % Colorbar
    hcb = colorbar;
    set(get(hcb,'ylabel'),'string',cb_label, 'FontSize', 12);
    if cb_auto == 0 
        caxis([0 32]);
    end

    % Title Etc..
    grid on;
    title_date = sprintf('%s - %s', datestr(start_date, 'mm/dd/yyyy'), ...
                 datestr(datenum(start_date+n_days-1), 'mm/dd/yyyy'));
    title_grid = sprintf('%s %d', names{i}, grid_num);
    title({title_grid;title_date}, 'FontSize', 12);

    % Anti-alias
    % aa = myaa(4);

    % Save image as *.png to post
    image_str = sprintf('%s_%s_%d_%s_%i.png', names{i}, model_variable, ...
                grid_num, datestr(start_date, 'mm-dd'), n_days);
    saveas(fH,image_str);
end