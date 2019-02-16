index = 1;
for station = 1:size(stations,2)
    for month = 1:12
        if ~isnan(tg(station).monthlyElev(month))
            elev(index,:) = [tg(station).lat, tg(station).monthlyElev(month),     ...
                             db14(station).lat, db14(station).monthlyElev(month), ...
                             db16(station).lat, db16(station).monthlyElev(month), ...
                             db22(station).lat, db22(station).monthlyElev(month), ...
                             month];
            index = index + 1;
        end
    end
end

