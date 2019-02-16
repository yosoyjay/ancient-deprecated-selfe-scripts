function [] = bndry_vert_profile(start_date, time_steps )
% This function plots a vertical profile of a node on the ocean open boundary.
%
% Input:
% start_date - The date that the run begins
% time_steps - The number of time steps to extract values
% model_variable - The model output variable to plot
%

v_station = [ 29869.38, 310854.91 ];

% Check path and availability of modext
ret_code = exist('modext.m','file');
if ret_code ~= 2
	addpath '/usr/local/cmop/modext';
end

% Data collected every 15 minutes at every 0.1 m to 25 meters
% n_days*output_dt - 1 because that is the actual number of available times.
% output_dt - The number of output time steps per day
output_dt = 96 
start_date = datenum(start_date);
times = [];
for i = 1:time_steps
	times = [times datenum(start_date + (1/output_dt)*i)];
end

depth = 2782;
depths = 0:(depth/100):depth;

model_output = modext(times, v_station(1,2), v_station(1,1), depths, '1_trcr_1.63', './');
imagescnan(model_output);

% Colorbar
hcb = colorbar;
set(get(hcb,'ylabel'),'string',cb_label, 'FontSize', 12);
if cb_auto == 0 
	caxis([0 32]);
end

% Title Etc...

% Save image as *.png to current directory
saveas(gcf, 'trcr_at_bndry.png');
