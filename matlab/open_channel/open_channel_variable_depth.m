function [] = open_channel_variable_depth(time)
% This function extracts data from open channel runs from the GLS Warner experiment. 
%
% Input:
% time - Time to extract data in datenum format
%
% Output:
% image(s) - Of variables as *.png files
% data     - Of variables as *.mat files
%
% lopezj 9/16/2011

% map file names to plot labels 
file_key = {'hvel.64', 'vdff.63', 'tdff.63', 'conc.63',   ...
            'kine.63', 'dahv.62'}; %, 'diss.63'};
label_val = {'Velocity, m/s', 'Eddy viscosity, m^2 s^{-1}', ...
			 'Eddy diffusivity, m^2 s^{-1}', 'Density, kg m^{-3}' ...
			 'Turbulence kinetic energy, m^2 s^{-2}', 'Depth Averaged Velocity, m/s'}; % ...
%			 'dissipation, m^2 s^-3'};
%labels = containers.Map(file_key, label_val);
%cb_label = labels(model_variable);

% Location to get data - 1km from boundary condition
v_station = [100, 500]; 


% Check path and availability of modext
ret_code = exist('modext.m','file');
if ret_code ~= 2
	addpath '/usr/local/cmop/modext';
end

% Data collected once during steady state every 0.5 m to 10 meters
%start_date = datenum(time);
%times = [];
%output_dt = 96;
%n_days = 2;
%for i = 1:(n_days*output_dt-1)
%	times = [times datenum(start_date + (1/output_dt)*i)];
%end
times = datenum(time);
depths = 0:0.5:10;
% Loop over every variable from the model and extract the data
for i = 1:size(file_key,2)
		model_results = modext(times, v_station(1,2), v_station(1,1), depths, file_key{i},'./');
		
		% Create plot 
		plot(model_results, depths);
		set(gca,'ydir','reverse');

		% X-Axis
		%set(gca,'xtick', 0:output_dt:(n_days*output_dt));
		%day_one = datestr(datenum(start_date+1),'mm/dd');
		%day_two = datestr(datenum(start_date+2), 'mm/dd');
		%xtick_labels = {day_two, day_one};
		%set(gca,'xticklabel', label_val{1});

		% Y-Axis
		%set(gca,'ytick', 0:20:101);
		%set(gca,'yticklabel', depths(1:20:101));
		ylabel('Depth From Free Surface (m)', 'FontSize', 12);

		% Title Etc..
		%grid on;
		%title_date = sprintf('%s - %s', datestr(start_date, 'mm/dd/yyyy'), ...
		%             datestr(datenum(start_date+n_days-1), 'mm/dd/yyyy'));
		%title_grid = sprintf('%s %s', names{i}, grid_num);
		title(label_val{i}, 'FontSize', 12);

		% Anti-alias
		% aa = myaa(4);

		% Save image as *.png to post
		image_str = sprintf('%s.png', file_key{i}); 
		saveas(gcf,image_str);

		% Save data from model_results
		save(file_key{i}, 'model_results')
end
