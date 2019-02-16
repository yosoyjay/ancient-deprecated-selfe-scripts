function [runInfo] = get_run_info(runPath)
% Returns infomation about a model run in the runPath
%
% Input:
%	runPath - Path to run directory
% Output:
%	runInfo - Structure with:
%		nSteps    - Number of steps in a run file
%		dt     	  - Model output timestep
%       startDate - Model start date in datenum
%		nDays	  - Number of days in model run
%
% Assumes that ?_elev.61 is output and works for my density project.
%
% lopezj - 1/24/2012
%

% Get information about start date
runInfo.startDate = get_start_date(runPath);

% Dates to datenum and interpolated to timestep of model
% Get dt and nSteps from '1_elev.61'
filePath = sprintf('%s/outputs/5_elev.61', runPath);
hElev = sz_readHeader(filePath);
runInfo.dt = hElev.dt;
runInfo.nSteps = hElev.nSteps;

% Number of days based on file number of last combined file 
% Weak sauce and good only for me with days < 10
cmd = sprintf('ls -1 %s/outputs/?_elev.61', runPath);
[s, r] = unix(cmd);
runInfo.nDays = str2num(r(end-9)); 

