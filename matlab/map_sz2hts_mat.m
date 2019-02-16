function [data_hts] = map_sz2hts_mat(eb_h,eb_data)
%function [data_hts] = map_sz2hts(eb_h, data_eb)
% map a <del>vector<\del> maxtrix of data as extracted from elcirc or selfe (sz) binary files dat(kbp:nvrt,1:np,1:ts)
% to array including filling in of the values above and below kfp kbp with nans.
%
% Input:
% 	eb_h    - A header for the binary file returned by ??_readHeader
% 	eb_data - The data read by ??_readTimeStep 
%
% Output:
%	data_hts - Data mapped to hts of the form (1:np,1:nvrt,1:ts) where ts is number of 
%			   output timesteps in eb_data
%
% Sergey Frolov March 2004
% modified from map_eb2hts_nan by SF, Jan 2006
% Modified to work with matrices by lopezj, 02/2012
%
data_hts = nan(size(eb_h.idx.idx_all), size(eb_data,3));
data_hts(eb_h.idx.idx_all_mask,:)=eb_data(eb_h.idx.idx_all(eb_h.idx.idx_all_mask),1,:);
data_hts = reshape(data_hts,eb_h.hgrid.np,eb_h.vgrid.nLevels, size(data_hts,2));
%data_hts(~eb_h.idx.idx_all_mask)=nan;
