function [row, col] = get_ncom_indices(lat, long, ncomLat, ncomLong)
% Returns the indices of NCOM data in NetCDF files closest to the lat and long
% given in the args. Uses Euclidean distance as the distances should be small
% enough to assume an approximate plane.
%
% I use it as part of a workflow to compare NCOM to observation data.
%
% Input:
%	lat      - Latitude to find
%	long     - Longitude to find
%   ncomLat  - Array of latitudes returned from get_ncom_ssh
%	ncomLong - Array of longitudes returned from get_ncom_ssh
%
% Output:
%   row 	 - Row indices into NCOM data
%   column   - Column indices into NCOM data
%
% lopezj - 01/03/2012
%

addpath '/home/workspace/users/lopezj/scripts/matlab';

% Flatten NCOM matrices
latLongSize = size(ncomLat,1)*size(ncomLong,2);
latFlat  = reshape(ncomLat,latLongSize,1);
longFlat = reshape(ncomLong,latLongSize,1);

% Use nearest neighbor to find closet match
idx = nearestneighbour([lat; long], [latFlat'; longFlat']);

% Return indices into the data's original shape 
[row, col] = ind2sub(size(ncomLat),idx);
