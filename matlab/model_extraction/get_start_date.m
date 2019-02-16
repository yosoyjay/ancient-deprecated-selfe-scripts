function [startDate] = get_start_date(runPath)
% Returns the start date of a run in datenum format
%
% Input:
%	runPath   - Path to run directory
% Output:
%	startDate - Start date in datenum format
%
% lopezj - 01/19/2012
%

% Get information about the start date from bctides - lifted from Grant's stuff
bctidesPath = sprintf('%s/bctides.in', runPath);
fid = fopen(bctidesPath, 'r');
dstr = fgets(fid);
startDate = datenum(dstr(1:end-4),'mm/dd/yyyy HH:MM:SS');
fclose(fid);

