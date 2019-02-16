function [ncomBias] = calc_ncom_bias(ctdData, ncomData, minDepth)
% This function calculates the bias of NCOM salinity data via the following
% method:
%  1. Calculate residual between ncom and smoothed ctd (ncom-ctd) span of 50
%	  - ctd data is interpolated (really kinda subsampled) to the depths of NCOM
%  2. Sum over the residuals to get a single residual per ctd cast
%
% Input:
%  ctdData - ctdData in the structure created by get_ctd_data
%  ncomData - ncomData in the structure created by get_ncom_data
%  minDepth - Minimum depth to consider, used to avoid freshwater effects in water column
% Output:
%  ncomBias[size(ctdCast,2)] - An array that contains the ncomBias for each cast.
%
% lopezj - 01/12/2012
%

ncomBias = nan(size(ctdData,2),1);

% Loop over every ctdData in the structure
for cast = 1:size(ctdData,2)
	% Find where data ends in NCOM or ctd data and interp ctd salinity to NCOM depths
	idx = find(isnan(ncomData(cast).data(:,1)),1)-1;	
	maxObsDepth = -max(ctdData(cast).data(:,1));
	if idx > 0
		if maxObsDepth > minDepth
			continue;
		end
		if maxObsDepth > ncomData(cast).data(idx,1)
			idx = find(ncomData(cast).data(:,1) < maxObsDepth,1);
		end
	else
		continue;
	end

	% Smooth ctd data
	sSmooth = smooth(ctdData(cast).data(:,2),50);
	dSmooth = smooth(ctdData(cast).data(:,1),50);

	% Interp ctd salinity to depths available in NCOM data
	sInterp = interp1(-dSmooth, sSmooth, ncomData(cast).data(1:idx,1), 'nearest','extrap'); 

	% Find residual at every NCOM depth and calculate the sum
	ncomBias(cast) = sum(ncomData(cast).data(1:idx,2)-sInterp);
end
