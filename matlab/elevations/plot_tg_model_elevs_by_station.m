year = 2000;
monthName = {'January', 'February', 'March', 'April', 'May', 'June', 'July', ...
             'August', 'September', 'October', 'November', 'December'};
stationName = {'South Beach', 'Garibaldi', 'Tongue Point', 'Toke Point', ...
               'Westport', 'La Push', 'Neah'};

% Loop over elvation data and derive indices into data to divide data by station 
plotIndexing = NaN(7,2);
i = 1;
j = 2;
while j <= size(elev,1)
    if elev(i,1) == elev(j,1)
        j = j+1;
        % Ugly, but I missed an edge case first time around
        if j == size(elev,1)
            switch elev(i,1)
                case latitudes(1)
                    x = 1;
                case latitudes(2)
                    x = 2;
                case latitudes(3)
                    x = 3;
                case latitudes(4)
                    x = 4;
                case latitudes(5)
                    x = 5;
                case latitudes(6)
                    x = 6;
                case latitudes(7)
                    x = 7;
            end
            plotIndexing(x, 1) = i;
            plotIndexing(x, 2) = j;
        end
    else
        switch elev(i,1)
            case latitudes(1)
                x = 1;
            case latitudes(2)
                x = 2;
            case latitudes(3)
                x = 3;
            case latitudes(4)
                x = 4;
            case latitudes(5)
                x = 5;
            case latitudes(6)
                x = 6;
            case latitudes(7)
                x = 7;
        end
        plotIndexing(x, 1) = i;
        plotIndexing(x, 2) = j-1;
        i = j;
        j = j+1;
    end
end

% Hacked up to get SSH from NCOM to work.
ssh_range = 1:12;
% Maps stations index to NCOM index specific to 4 stations for 2000
ssh_index = [1 NaN 2 3 NaN NaN 4];  
% Loop over stations, if plotIndexing has values plot it 
for i = 1:7
    clf;
    if ~isnan(plotIndexing(i,1));
        range = plotIndexing(i,1):plotIndexing(i,2); 
        plot(elev(range,9), elev(range,2), 'color', 'r', 'marker', 'o', 'linewidth', 1); hold on;
        plot(1:12, db14(1,i).monthlyElev(1,:), 'color', 'b', 'marker', 'x', 'linewidth', 1); hold on;
        plot(1:12, db16(1,i).monthlyElev(1,:), 'color', 'g', 'marker', '+', 'linewidth', 1); hold on;
        plot(1:12, db22(1,i).monthlyElev(1,:), 'color', 'y', 'marker', 's', 'linewidth', 1); hold on;
        plot(ssh_range, meanSSH(ssh_index(i),:), 'color', 'k', 'marker', '*', 'linewidth', 1); hold on;
%       plot(elev(range,5), elev(range,2)-elev(range,4), 'color', 'k'); hold on;
        legend('Tide Gauge', 'DB14', 'DB16', 'DB22', 'NCOM');
        ylabel('Elevation m');
        ylim([-0.3 0.6]);
        xlabel('Month');
        xlim([1 12]);
        title(stationName{i});
        figureName = sprintf('%02d_%04d_station.png', i, year); 
        saveas(gcf, figureName);                
        figureName = sprintf('%02d_%04d_station.fig', i, year); 
        saveas(gcf, figureName);                
        clf;
    end
end
