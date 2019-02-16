classdef modelStation
    % Model station class
    properties
        name        % Station name
        lat         % Latitude (decimal)
        long        % Longitude (decimal)
        number      % Station number. 1 being to most Southern, n being the most Northern
        year        % Year of this data
        db          % DB of the model (11, 14, 22, etc.)
        dataPath    % Path to directory of data
        rawData     % Raw data from boundary nodes
        monthlyElev % Monthly averaged Elevation
    end

    methods
        function obj = modelStation(number, name, latitude, longitude, db, year);
            obj.name      = name;
            obj.number    = number;
            obj.lat       = latitude;
            obj.long      = longitude;
            obj.db        = db;
            obj.year      = year;
            obj.monthlyElev = NaN(1,12);
        end
        function obj = calcElev(obj, month)
            dateString = sprintf('%d-%02d-01 01:00:00', obj.year, month);
            startDate = datenum(dateString);
            date(1) = startDate;
            for time=2:(24*30);
                date(time) = addtodate(date(time-1), 1, 'hour');
            end
            db = sprintf('DB%d', obj.db);
            obj.monthlyElev(month) = mean(modext(date, obj.lat, obj.long, 0, 'Elev', db));          
        end 
    end 
end


