function [ctdData] = get_ctd_data(filePath, saveAsMat)
% Returns a struct with CTD data or an array of structs
% if there is more than one cast per file
%
% Inputs:
%   filePath - Path to the *.mat ctd data
%   saveAsMat - Binary, 1 - save as *.mat, 0 - do not save 
% Outputs:
%   ctdData - Struct of CTD data
%       lat
%       long
%       date
%       data[] - depth, salinity, temperature
%
% lopezj - 01/03/2010
%

% Load all of the data from the CTD database file and first line
data = importdata(filePath,',');
tmpCastID = data.textdata(2,1);
ctdIdx = 1;
ctdData(ctdIdx).castID = tmpCastID;
ctdData(ctdIdx).lat    = data.data(1,3);
ctdData(ctdIdx).long   = data.data(1,4);
tempTime = data.textdata{2,3}(1:(end-3));  % The only reason this is here is to prevent Matlab from chocking on the next line.
ctdData(ctdIdx).time   = datenum(tempTime);
ctdData(ctdIdx).data   = [str2num(cell2mat(data.textdata(2,2))) data.data(1,1) data.data(1,2)];

% Seperate data into diff structs for each cast ID
for i = 2:size(data.data,1)
    if ~strcmp(cell2mat(tmpCastID), cell2mat(data.textdata(i+1,1)))
        ctdIdx = ctdIdx + 1;
        tmpCastID = data.textdata(i+1,1);
        ctdData(ctdIdx).castID = tmpCastID;
        ctdData(ctdIdx).lat    = data.data(i,3);
        ctdData(ctdIdx).long   = data.data(i,4);
        ctdData(ctdIdx).time   = datenum(data.textdata{i+1,3}(1:(end-3)));
        ctdData(ctdIdx).data   = [str2num(cell2mat(data.textdata(i+1,2)))...
                                  data.data(i,1) data.data(i,2)];
    else
        ctdData(ctdIdx).data   = [ctdData(ctdIdx).data; ...
                                  str2num(cell2mat(data.textdata(i+1,2)))...
                                  data.data(i,1) data.data(i,2)];
    end
end

% Save *.mat if specified 
if saveAsMat 
    matFile = sprintf('ctdData-%s.mat', datestr(ctdData(1).time, 'yyyy'));
    save(matFile, 'ctdData');
end
