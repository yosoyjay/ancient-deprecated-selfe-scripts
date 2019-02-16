function [] = vertical_profile_animation(start_date, n_days, model_variable, grid_num, node_file)
%
% This function makes a vertical profile from model results at saturn1
%
% Input:
% start_date - First day to extract data 
% n_days     - Number of days to extract
% model_variable - What variable from the model do you want? 'Salt', 'elevation', etc. 
% grid_num   - What grid is used? (Only for plot extraction)
% node_file  - File that has the node numbers to extract that data from
%
% Output: Nones
%
% lopezj - 6/19/2011
%

% Number of times to get data per day. 96 is Max.
    DAY_DIV = 48;

% Check path and availability of modext
    ret_code = exist('modext.m','file');
    if(ret_code ~= 2)
        addpath '/usr/local/cmop/modext';
    end

% Get the positions of the nodes where you want data
    v_stations = get_node_cords(node_file, grid_num);

% n_days*DAY_DIV - 1 because that is the actual number of available times.
    start_date = datenum(start_date);
    times = [];
    for i = 1:(n_days*DAY_DIV-1)
        times = [times datenum(start_date + (1/DAY_DIV)*i)];
    end

    % Depths adequate to resolve deepest part of N. Channel
    depths = 0:0.25:25;
% Create matrices for data and plotting
    salinity = zeros(size(depths,2), size(v_stations,1));
    elevations = zeros(size(v_stations,1),1);
    corrected_depths = zeros(size(depths,2), size(v_stations,1));

% X-axis - A point for each station
    x = ones(size(depths,2), size(v_stations,1));
    for i = 1:size(v_stations,1)
        x(:,i) = x(:,i)*i;
    end
    v_stations_label = round(v_stations/1000);

% Loop over all times, create plot, and save image
    for i = 1:size(times,2)

        % Extract data you are looking for
        for j = 1:size(v_stations,1)
            %salinity = modext(times, 46.235, -123.869, depths, model_variable, './');
            salinity(:,j) = modext(times(i), v_stations(j,2), v_stations(j,1), depths, model_variable, './');
        end

        % Get elevation to corrrect depth
        for j = 1:size(v_stations,1)
            elevations(j,1) = modext(times(i), v_stations(j,2), v_stations(j,1), 0, 'Elev', './');
        end

        % Correct depth
        for j = 1:size(v_stations,1)
            corrected_depths(:,j) = depths - elevations(j);
        end

        % Create image and labels
        %imagesc(salinity);
        %set(gca,'xtick', 0:1:size(v_stations,1));
        if i == 1
            h = gcf;
            a = imgca(h);

            % General Plot Format - Highlight Saturn01 Position
            view(a, 0,-90);
            annotation(h,'textbox',...
                [0.566202090592335 0.924468493397226 0.113240418118467 0.0656565656565656],...
                'String',{'Saturn 1'},...
                'FitBoxToText','off',...
                'LineStyle','none');
            annotation(h,'line',[0.608429484063389 0.608429484063388],...
                [0.108287777552811 0.928994848259882],'LineWidth',2);

            % Y-Axis
            ylabel(a, 'Depth (m)');
            set(a,'ytick', -5:5:25);
            set(a,'yticklabel', -5:5:25);

            % X-Axis 
            set(a,'xtick', 1:5:40);
            set(a,'xticklabel', v_stations_label(1:5:40));
            xlabel(a,'Longitude-SPCS (m)');

            % Colorbar
            hcb = colorbar;
            set(get(hcb,'ylabel'),'string','Model Salnity (psu)');
            caxis([0 32]);

            % Title
            title_date = sprintf('%s', datestr(times(i), 'mm/dd/yyyy HH:MM'));
            title_grid = sprintf('%s %i', 'North Channel', grid_num);
            title(a, {title_grid;title_date});

        end

        % Generate the plot
        p = surf(a, x, corrected_depths, salinity, 'EdgeColor', 'none');
    
        % Save image as *.png to post
        image_str = sprintf('%03d_%s_%i_%s_%i.png', i, model_variable, grid_num, datestr(start_date, 'mm-dd'), n_days);
        saveas(h,image_str);

        % Clean slate 
        delete(p);
    end
end
