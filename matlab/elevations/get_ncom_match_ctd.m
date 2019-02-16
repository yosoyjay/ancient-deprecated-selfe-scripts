function [ncomData] = get_ncom_match_ctd(ctdData)
% Given a set of ctdData created by get_ctd_data, this script will find NCOM
% data that is closest spatially and temporally and calculate the residual.
%
% I used this to compare NCOM to observation data
%
%   
% Input:
%   ctdData - Array of ctdData structs created by get_ctd_data
% Output:
%   ncomData - Array of NCOM structs closest to observed data
%
% lopezj - 01/04/2012
%

% All levels in NCOM model
depthRange = 1:39;

% For every ctd cast...
for i = 1:size(ctdData,2)
    % Get NCOM data for the day of the CTD cast
    [salt, temp, depth, ncomLat, ncomLong] = get_ncom_salTemp(datevec(ctdData(i).time));
	salt(salt < 0) = NaN;
	temp(temp < 0) = NaN;
	depth(depth < -5000) = NaN;
    % Get the indices for the closest NCOM node to the CTD cast
    [row, column] = get_ncom_indices(ctdData(i).lat, ctdData(i).long, ncomLat, ncomLong);
    % Create ncomData structure
    ncomData(i).lat  = ncomLat(row, column);
    ncomData(i).long = ncomLong(row, column); 
    ncomData(i).time = ctdData(i).time;
	ncomdData(i).data = zeros(39,3);
	for j = 1:depthRange(end)
    	ncomData(i).data(j,1) = depth(row, column, j); 
		ncomData(i).data(j,2) = salt(row, column, j); 
		ncomData(i).data(j,3) = temp(row, column, j);
	end
	ncomData(i).row = row;
	ncomData.column = column;
end

    
