function [assessment,lines,mod_data,obs_data] = skillassessment2(varargin)
% [assessment,lines,mod_data,obs_data] = SKILLASSESSMENT2(run_type,mod_name,year,[doy day] [week day],outfile,basedir)
%
% Calculates skill assessment indices for the given forecast [hindcast]
% using the skill assessment m-files currently in the /safunc_modules
% file.
% 
% run_type = 'forecast' or 'hindcast', which is used to interpret arguments
% mod_name = string name of the simulation as it appears in the
%                         directory call (i.e. 'dev', 'db16', '14')
% year = the year in which the simulation to assess occurs
% week = the week in which the simulation to assess occurs (used
%               only for hindcasts)
% day = the day of year of the forecast
% day (optional) = the day of the run (1-7 for hindcasts, 
%             1-2 or 1-3 for forecsasts), if omitted, the entire run is processed
% outfile (optional) = location to write skillassessment data to, argument of '' means don't write file, argument of 0 means use default name (specify if specifying basedir)
% basedir (optional) = location of run directories (include full path below the year-week-mod_name or year-doy level, including mod_name in forecasts), defaults to '/home/workspace/ccalmr/hindcasts/' for hindcasts and to ['/home/workspace/local0/forecasts/' mod_name] for forecasts
% example: hcast_assess = skillassessment('hindcast','14',2006,[24 5]); skill assessment for the 5th day of hindcast 2006-24-14
%          hcast_assess = skillassessment('hindcast','14',2006,24); skill assessment for the full week of hindcast 2006-24-14
%          fcast_assess = skillassessment('forecast','db16',2006,182,0,'home/workspace/ccalmr/forecasts/dev/2008/'); for a forecast, using a non-standard directory, with the standard output file name
% lines = contents of assessment rearranged for printing to a file
sim = varargin{1};
mod_name = varargin{2};
year = varargin{3};
weekday = varargin{4};
if length(weekday) == 2
   day = weekday(2);
else
   day = 0;
end
if nargin == 6
   basedir = varargin{6};
else
   basedir = '';
end
if nargin >=5
   outfile = varargin{5};
else 
   outfile = 0;
end
% set up forecast/hindcast specific variables and gather data
switch 1
    case strncmp(sim,'f',1)
        sim = 'forecast';
        doy = weekday(1);
        if isempty(basedir)
            basedir = ['/home/workspace/local0/forecasts/' mod_name ,'/'];
        end
        rundir = [basedir, num2str(year) '-' sprintf('%03d',doy) '/'];
% gather the appropriate data
        [mod_data,obs_data,goodsta] = process_forecast(mod_name,year,doy,day,basedir);
    case strncmp(sim,'h',1)
        sim = 'hindcast';
        week = weekday(1);
        if isempty(basedir)
            basedir = '/home/workspace/ccalmr/hindcasts/';
        end
        rundir = [basedir, num2str(year) '-' sprintf('%02d',week) '-' mod_name '/'];
% gather the appropriate data
        [mod_data,obs_data,goodsta] = process_hindcast(mod_name,year,week,day,basedir);
    case varargin<4 | varargin>6
        help skillassessment
end
if outfile == 0   % if outfile not specified at command line, generate default outfile
  if day > 0  % outfile includes day of model, unless the skill assessment is done onthe full run
     outfile = [rundir 'process/skillassessment_' sprintf('%03d',day) '.csv'];
  else 
     outfile = [rundir 'process/skillassessment.csv'];
  end
end

% Filter datasets
%   -Remove zero salinity pairs
for sta = 1:goodsta
    for var = 1:length(mod_data(sta).variable)
        if strcmp(mod_data(sta).variable(var).name(1:4),'sali')
            bad = find(mod_data(sta).variable(var).data==0&obs_data(sta).variable(var).data==0);
            mod_data(sta).variable(var).data(bad) = [];
            obs_data(sta).variable(var).data(bad) = [];
        end
    end
end

% Determine which index modules are available
%funclist =  dir('safunc_modules/*calc.m');
%funclist =  dir('/usr/local/matlab/toolbox/cmop/skill_assessment/safunc_modules/*calc.m');
funclist = dir('/home/choj/matlab/skill_assessment/trunk/safunc_modules/*calc.m');
% addpath safunc_modules/;
for i = 1:length(funclist)
    m_file = strtrim(funclist(i).name);
    module{i} = m_file(1:end-2);
end
modules = length(module);
% Calculate indices
% First do calculations at broadest scale
%model = [];
%observation = [];
%for i = 1:length(mod_data)
%    for j = 1:length(mod_data(i).variable)
%        model = [model;mod_data(i).variable(j).data];
%        observation = [observation;obs_data(i).variable(j).data];
%    end
%end
if length(mod_data) > 0
%for ind = 1:length(module)
%    eval(['assessment.' module{ind}(1:end-4) ' = ' module{ind} '(observation,model);'])
%end
%assessment.n = length(model);
% Now do at the station scale
%for i = 1:length(mod_data)
%    model = [];
%    observation = [];
%    for j = 1:length(mod_data(i).variable)
%        model = [model;mod_data(i).variable(j).data];
%        observation = [observation;obs_data(i).variable(j).data];
%    end
%    for ind = 1:modules
%        eval(['assessment.station(i).' module{ind}(1:end-4) ' = ' module{ind} '(observation,model);'])
%        assessment.station(i).name = mod_data(i).name;
%        assessment.station(i).n = length(model);
%    end
%end
% Now do at the station/variable scale
var_name = {};
cnt = 0;
assessment = struct('n',0,'station',struct('n',0,'variable',struct),'variable',struct('n',0),'weight',0);
for i = 1:length(mod_data)
    assessment.station(i).n = 0;
    assessment.station(i).weight = 0;
    for j = 1:length(mod_data(i).variable)
        new = 1;
        for k = 1:length(var_name)
            if strcmp(var_name{k}(1:4),mod_data(i).variable(j).name(1:4))
                new = 0;
                varind = k;
            end
        end
        if new==1
            cnt = cnt+1;
            varind = cnt;
            var_name{cnt} = mod_data(i).variable(j).name;
            assessment.variable(varind).n = 0;
            assessment.variable(varind).name = mod_data(i).variable(j).name;
            assessment.variable(varind).weight =  0;
        end
        model = mod_data(i).variable(j).data;
        observation = obs_data(i).variable(j).data;
        assessment.station(i).variable(j).n = 0;
        assessment.station(i).n = assessment.station(i).n + length(model);
        try 
           assessment.variable(varind).n = assessment.variable(varind).n + length(model);
        catch
           keyboard
        end
        assessment.n = assessment.n + length(model);
        
        assessment.station(i).variable(j).name = mod_data(i).variable(j).name;
        assessment.station(i).variable(j).n = length(model);
        w = length(model).*mod_data(i).dt;
        assessment.station(i).variable(j).weight = w;
        assessment.station(i).weight = assessment.station(i).weight + w;
        assessment.variable(varind).weight =  assessment.variable(varind).weight + w;
        assessment.weight =  assessment.weight + w;
        for ind = 1:modules
             modname = module{ind}(1:end-4);
            eval(['assessment.station(i).variable(j).' module{ind}(1:end-4) ' = ' module{ind} '(observation,model);'])
            if isfield(assessment.station(i),modname)
            	if eval(['~isempty(assessment.station(i).' modname ');']);
            	   eval(['assessment.station(i).' modname ' = assessment.station(i).' modname ' + assessment.station(i).variable(j).' modname ' * w;']);
	        else 
                   assessment.station(i).name = mod_data(i).name;
                   eval(['assessment.station(i).' modname ' = assessment.station(i).variable(j).' modname ' * w;']);
                end
	    else 
                assessment.station(i).name = mod_data(i).name;
                eval(['assessment.station(i).' modname ' = assessment.station(i).variable(j).' modname ' * w;']);
            end
            if isfield(assessment.variable(varind),modname)
                if eval(['~isempty(assessment.variable(varind).' modname ');'])
            	    eval(['assessment.variable(varind).' modname ' = assessment.variable(varind).' modname ' + assessment.station(i).variable(j).' modname ' * w;']);
	         else 
                    eval(['assessment.variable(varind).' modname ' = assessment.station(i).variable(j).' modname ' * w;']);
                 end
	    else 
                eval(['assessment.variable(varind).' modname ' = assessment.station(i).variable(j).' modname ' * w;']);
            end
            if isfield(assessment,modname)
                if  eval(['~isempty(assessment.' modname ');'])
            	    eval(['assessment.' modname ' = assessment.' modname ' + assessment.station(i).variable(j).' modname ' * w;']);
	        else 
                    eval(['assessment.' modname ' = assessment.station(i).variable(j).' modname ' * w;']);
                end
	    else 
                eval(['assessment.' modname ' = assessment.station(i).variable(j).' modname ' * w;']);
            end
        end
        %keyboard
    end
end
ass =assessment;
for k = 1:modules
   modname = module{k}(1:end-4);
   eval(['assessment.' modname  ' =  assessment.' modname  ' ./ assessment.weight;']); 
   for i=1:length(assessment.station)
       eval(['assessment.station(i).' modname  ' =  assessment.station(i).' modname  ' ./ assessment.station(i).weight;']);
   end
   for i=1:length(assessment.variable)
      eval(['assessment.variable(i).' modname  ' =  assessment.variable(i).' modname  ' ./ assessment.variable(i).weight;']);
   end
end
% Calculate perfect scores for the current set of indices
obs = ones(10,1);
mod = ones(10,1);
for ind = 1:modules
    eval(['assessment.perfect(ind).value = ' module{ind} '(obs,mod);'])
    assessment.perfect(ind).name = module{ind};
end
lines = sa_preplines(assessment);
outfile
if ~isempty(outfile)
  f=fopen(outfile,'w');
  for i=1:length(lines)
    fprintf(f,'%s\n',lines{i});
  end
  fclose(f);
end
else 
assessment = struct();
lines = {};
end
%------------------------------------

function [mod_data,obs_data,goodsta] = process_forecast(mod_name,yr,dy,rday,basedir)
% Gather observational data;
year = num2str(yr);
day = num2str(dy,'%03d');
if rday > 0
   rday = num2str(rday,'_n%03d');
else
   rday = '';
end
data_stations = 0;

% Find appropriate files
list = dir([basedir, '/' year '-' day '/data/*CTD']);
station = {list.name};
if isempty(station)
    error(['Couldn''t find any station data in ' basedir, '/' year '-' day '/data/*CTD'])
end
stations = length(station);
%ctdfiles = [];
%keyboard
%for sta = 1:stations
%    ctdfiles(sta) = strcmp(station{sta}(end-2:end),'CTD');
%end
%stations = sum(ctdfiles);
%station_old = station;
%station(~ctdfiles) = [];

for sta = 1:stations
    underscore = findstr('_',station{sta});
    statname{sta} = strtrim(station{sta}(1:underscore(1)-1));
    instr{sta} = strtrim(station{sta}((underscore(1)+1):(end-8)));
end

% Create list of filenames for model data
goodsta = 0;
mod_data = [];
obs_data = [];
for sta = 1:stations
    data_stations = data_stations+1;
    dstatname{data_stations} = strtrim(statname{sta});
    statname{sta}
    depth = sa_getdepth(statname{sta},instr{sta})
	mod_day_dir = [basedir, '/'  year '-' num2str(day) '/process/' statname{sta} '_elcirc' rday '.dat'];
	obs_day_dir = [basedir, '/'  year '-' num2str(day) '/data/' station{sta}];
        [mod_new,obs_new,dt] = sa_getdata(mod_day_dir,obs_day_dir,depth);
    % Store data in the data structures
	if ~isempty(mod_new)
        goodsta = goodsta+1;
        mod_data(goodsta).variable = [];
        mod_data(goodsta).dt = dt;
        mod_data(goodsta).name = [dstatname{data_stations} '_' instr{sta} '_' sprintf('%05d',-depth*100)];
        obs_data(goodsta).variable = [];
        obs_data(goodsta).dt = dt;
        obs_data(goodsta).name = [dstatname{data_stations} '_' instr{sta} '_' sprintf('%05d',-depth*100)];
        for datatype = 1:length(mod_new)
            mod_data(goodsta).variable(datatype).name = mod_new(datatype).name;
            mod_data(goodsta).variable(datatype).data = mod_new(datatype).data;
            mod_data(goodsta).variable(datatype).t = mod_new(datatype).t;
            obs_data(goodsta).variable(datatype).name = obs_new(datatype).name;
            obs_data(goodsta).variable(datatype).data = obs_new(datatype).data;
            obs_data(goodsta).variable(datatype).t = obs_new(datatype).t;
        end
	end
end

%------------------------------------

function [mod_data,obs_data,goodsta] = process_hindcast(mod_name,yr,wk,dy,basedir)

% Gather observational data;
year = num2str(yr);
week = sprintf('%02d',wk);
day = num2str(dy);
data_stations = 0;

% Hindcast data is stored in a central repository
databasedir = '/home/workspace/ccalmr/data/verified/model_weeks/';

% Find appropriate files
list = dir([databasedir, '*_' year '_w' week '.dat_*']);
station = {list.name};
if isempty(station)
    error('Couldn''t find any station data!')
end

ctdfiles = [];
for sta = 1:length(station)
    ctdfiles(sta) = strcmp(station{sta}(end-2:end),'CTD');
end
stations = sum(ctdfiles);
station_old = station;
station(~ctdfiles) = [];

for sta = 1:length(station)
    statname{sta} = strtrim(station{sta}(3:7));
    
end

if dy == 0
  dayfix = '';
else
  dayfix = sprintf('_n%03d',dy);
end
% Create list of filenames for model data
goodsta = 0;
mod_data = struct;
obs_data = struct;
for sta = 1:length(station)
    data_stations = data_stations+1;
    disp(station{sta})
    dstatname{data_stations} = strtrim(statname{sta});
    depth = -1*str2num(station{sta}(14:18))/100;
    instr = station{sta}(8:13);
	mod_day_dir = [basedir,num2str(year) '-' num2str(week) '-' mod_name '/process/' statname{sta} '_elcirc' dayfix '.dat'];
	obs_day_dir = [databasedir,station{sta}];
	[mod_new,obs_new,dt] = sa_getdata(mod_day_dir,obs_day_dir,depth);
    % Store data in the data structures
	if ~isempty(mod_new)
        goodsta = goodsta+1;
        mod_data(goodsta).variable = [];
        mod_data(goodsta).dt = dt;
        mod_data(goodsta).name = [dstatname{data_stations} '_' instr '_' sprintf('%05d',-depth*100)];
        obs_data(goodsta).variable = [];
        obs_data(goodsta).dt = dt;
        obs_data(goodsta).name = [dstatname{data_stations} '_' instr '_' sprintf('%05d',-depth*100)];
        for datatype = 1:length(mod_new)
            mod_data(goodsta).variable(datatype).name = mod_new(datatype).name;
            mod_data(goodsta).variable(datatype).data = mod_new(datatype).data;
            mod_data(goodsta).variable(datatype).t = mod_new(datatype).t;
            obs_data(goodsta).variable(datatype).name = obs_new(datatype).name;
            obs_data(goodsta).variable(datatype).data = obs_new(datatype).data;
            obs_data(goodsta).variable(datatype).t = obs_new(datatype).t;
        end
	end
end
