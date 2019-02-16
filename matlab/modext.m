function [out1,out2] = modext(date,ypos,xpos,zpos,var,db_name,output_ts)
%MODEXT Model data extraction utility.
%   DATA = MODEXT(time,y,x,depth,variable,source), using the requested
%   database source, extracts the named variable at a particular
%   date, time, and depth. When the time, lat, lon, and depth arrays are
%   equal in size, MODEXT returns an array containing the variable values
%   for each (x,y,z,t). Up to two of the xyzt inputs (time, position, or
%   depth) may be passed as matrices. This is useful when creating time 
%   series, depth profiles, or transects. If no data are available at a
%   requested point and time, MODEXT will return NaN. Entering the string
%   'All' for time, x & y position, or depth, will return all values
%   available in the targeted binary data file.
%  
%   [OUT1,OUT2] = MODEXT(time,y,x,depth,4D_variable,source) passes the x
%   and y components of vector variables, such as velocity.
%
%   TIME may be expressed as a julian date or in string form
%   X,Y may be expressed in meters, or in latitude/longitude (dec. deg.)
%   DEPTH is expressed in meters below the free surface, 'Bot' for bottom
%         values, 'Surf' for surface values, or 'Ave' for the vertical
%         average
%   VARIABLE is either 'Salt', 'Temp', 'Elev', or 'Vel'. If you wish to
%         extract another variable, you can, but you must carefully enter
%         the variable name as it appears in the output files, including
%         the extension. For example, 'wind.62', or 'trcr_1.63'.
%   SOURCE is either the database to be extracted from (i.e., '14', '16',
%         or '22'), or the directory where a run's "outputs" folder is 
%         located. This option will only function accurately if the run's
%         bctides.in file is present, and contains a header line with the
%         correct date, expressed in the "mm/dd/yyyy hh:mm:ss PST" format.
%         No flexibility here, sorry.
%
%   EXAMPLE:
%      salinity = modext(733912.34,345000,290000,2.3,'Salt','14');
%      velocity = modext('2001-06-21 11:23:56',46.28,123.8,0,'Vel','16');
%      temperature = modext(733912.34,ones(1,11).*46.25, ...
%                           (124:.01:124.1),-2,'Temp','22');
%      [u,v] = modext(733912.34,46.28,123.8,(0:0.01:1),'Velocity','DB22');
%      elevation = modext('06/21/2001 11:23:56',46.28,123.8,0, ...
%                         'Elevation','16');
%      dates = ['2001-06-21 12:15:00'
%               '2001-06-21 12:30:00'
%               '2001-06-21 12:45:00'
%               '2001-06-21 13:00:00'];
%      salinity = modext(dates,46.28,123.8,0,'salinity','DB16');
%      tracer = modext(datenum(2002,3,22),46.28,123.8,0,'trcr_1.63', ...
%                      '/home/workspace/project/2010-07-14/run/');
%      salinity = modext(733912.463,'All','All','Surf','Salt','16');
%

%-------------------------------
%   Additional details:
%
%
%   MODEXT uses the compiled fortran codes developed by Joseph Zhang, and
%   the melio package developed by Sergey Frolov to directly access the
%   SELFE binary data files

%   C. Grant Law 09-10-2010

% First evaluate all the inputs and convert them to a standard format
%**************

% Look for remnants of past unsuccessful modext calls
if exist('modext_input.dat','file')
    !rm modext_input.dat
end
if exist('modext_output1.dat','file')
    !rm modext_output1.dat
end
if exist('modext_output2.dat','file')
    !rm modext_output2.dat
end

% Clean up and convert dates to julian days
jd = date_process(date);

% Determine if extraction will be M-Elio based
if  strcmp(date(1),'A')| ... % Conditions that indicate the use of melioext
    strcmp(date,'a') | ...
    isstr(xpos)| ...
    isstr(ypos)
    if isempty(strfind(db_name,'.6'))% If there are no dots, we need to find the location of the file
        [db,dir] = db_process(db_name);
        var = var_process(var);
        jd = date_process(date);
        zpos = depth_process(zpos);
        s = rem(jd,1).*86400;s(s==0)=86400;
        dv = datevec(jd); year=dv(:,1); yr = num2str(year); sol_day = floor(jd-datenum(year,1,1)+1);
        week = floor((sol_day-1)/7)+1; wk = num2str(week,'%.2i');
        iter = mod(sol_day-1,7)+1;
        for day = 1:length(jd)
            file = [num2str(iter(day)) '_' var];
            if ~isstr(db) % This file is from a calibration run
                runday = jd(day)-db;
                file = [num2str(ceil(runday)) '_' var];
                dbname(day,:) = [dir '/outputs/' file];
            else % This file is from a database run
                s = rem(jd,1).*86400;s(s==0)=86400;
                dv = datevec(jd); year=dv(:,1); yr = num2str(year); sol_day = floor(jd-datenum(year,1,1)+1);
                week = floor((sol_day-1)/7)+1; wk = num2str(week,'%.2i');
                iter = mod(sol_day-1,7)+1;
                file = [num2str(iter(day)) '_' var];
                dbname(day,:) = [dir yr(day,:) '-' wk(day,:) '-' db '/run/' file];
            end
            if ~exist(dbname(day,:),'file')
                error('Model results unavailable for that database/date')
            end
        end
    else % The target file has been specified
        if max(jd)-min(jd) > 1
            disp('When extracting from a single file, you')
            disp('may only specify dates within that day.')
        end
        dbname = db_name;
    end
    dind = false;
    first = 1;
    while ~dind(end)
        for i = 1:length(jd)
            dind(i) = strcmp(dbname(i,:),dbname(first,:));
        end
        in = find(dind); this = in(1); first = in(end)+1;
        out = fullfile_extract(dbname(this,:),date(dind),zpos);
        out1.x = out.x; out1.y = out.y; out1.bathy = out.bathy; out1.e = out.e;
        out1.data(dind) = out.data;
        out1.data(dind) = out.data;
        out1.jd(dind) = out.jd;
    end
    return
end

% Clean up and convert horizontal positions to model format
[x,y] = pos_process(ypos,xpos);

% Clean up and convert depth information
z = depth_process(zpos);

% Write input data to a common dimensional format
[jd,x,y,z] = redimensionalize(jd,x,y,z);

% Establish which variable is needed
varname = var_process(var);
dim = 1;
out1 = zeros(length(jd),1);
if strcmp(varname(end),'4')
    dim = 2;
    if nargout>1
        out2=out1;
    end
end

% Identify and clean up the name of the database to pull from
[db,dir] = db_process(db_name);

% Determine which code is needed, and modify arrays appropriately
[mode,jd,x,y,z] = mode_process(jd,x,y,z);

% Save matrix shape and rewrite input to linear arrays
block = zeros(size(jd));
jd = jd(:); x = x(:); y = y(:); z = z(:);

% Determine how many days are needed, and identify distinct dates
run_days = unique(floor(jd));
for day = 1:length(run_days)
    ind(:,day) = jd>=run_days(day)+output_ts/86400 & jd<=run_days(day)+1;
    a_ind(:,day) = jd>=run_days(day) & jd<run_days(day)+output_ts/86400;
end


%*************

% Now run the code
% Iterate through each of the days included in the observations
for day = 1:length(run_days)
    if any(ind(:,day)) % Are there any obs's in this day?
%     waitbar(day/(t1-floor(t0)),h);
        t = jd(ind(:,day));
        xset = x(ind(:,day));
        yset = y(ind(:,day));
        zset = z(ind(:,day));
        % Build input file
        build_input(t,xset,yset,zset,varname,db,dir,mode);
        % Get values
        eval(['!/usr/local/cmop/modext/' mode ' > /dev/null 2>&1'])
        !rm modext_input.dat
        % Assign to output
        data = load('modext_output1.dat');
        if any(data(:,2)==-99)
            data(data(:,2)==-99,:)=NaN;
        end
        out1(ind(:,day)) = data(:,2);
        !rm modext_output1.dat
        if dim==2
            data = load('modext_output2.dat');
            if any(data(:,2)==-99)
                data(data(:,2)==-99,:)=NaN;
            end
            if nargout==1
              out1(ind(:,day)) = sqrt(out1(ind(:,day)).^2 + data(:,2).^2);
            else
              out2(ind(:,day)) = data(:,2);
            end
            !rm modext_output2.dat
        end
  	end
    if any(a_ind(:,day)) % Are there any obs's in the first 15 minutes?
        % Build first input file
        t = jd(a_ind(:,day));
        xset = x(a_ind(:,day));
        yset = y(a_ind(:,day));
        zset = z(a_ind(:,day));
        s = rem(t,1).*86400;
        n = length(t);
        interp_data = zeros(n,1);
        val = [];
        % Determine which db22 files to access
        tfill = [floor(t(1))+output_ts/86400 floor(t(1))-.00001];
        for f=1:2
            % Build input file
            build_input(ones(size(xset)).*tfill(f),xset,yset,zset,varname,db,dir,mode);
            % Get values
            eval(['!/usr/local/cmop/modext/./' mode ' > /dev/null 2>&1'])
            !rm modext_input.dat
            % Assign to output
            data = load('modext_output1.dat');
            if any(data(:,2)==-99)
                data(data(:,2)==-99,:)=NaN;
            end
            val(:,f,1) = data(:,2);
            !rm modext_output1.dat
            if dim==2
                data = load('modext_output2.dat');
                if any(data(:,2)==-99)
                    data(data(:,2)==-99,:)=NaN;
                end
                val(:,f,2) = data(:,2);
                !rm modext_output2.dat
            end
        end
        % Now do the interpolation
        for it=1:dim
            interp_data = val(:,2,it)+s./output_ts.*(val(:,1,it)-val(:,2,it));
            if it==2
                if nargout==1
                    out1(a_ind(:,day)) = sqrt(out1(a_ind(:,day)).^2 + interp_data.^2);
                else
                    out2(a_ind(:,day)) = interp_data;
                end
            else
            	out1(a_ind(:,day)) = interp_data;
            end
        end
    end
end
%close(h);
block(:) = out1;
out1 = block;
if exist('out2','var')
    block(:) = out2;
    out2 = block;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function jd = date_process(date)

% Convert to julian days
if size(date,2)>1 && ~ischar(date)
    date = date';
end
switch true
    case ischar(date(1))
        switch true
            case length(findstr(date(1,:),'/'))==2 & any(findstr(date(1,:),'/')==5)
            	formstr = 'yyyy/mm/dd HH:MM:SS';
         	case length(findstr(date(1,:),'/'))==2 & any(findstr(date(1,:),'/')==3) & str2double(date(1,1:2))<13
            	formstr = 'mm/dd/yyyy HH:MM:SS';
          	case length(findstr(date(1,:),'-'))==2 & any(findstr(date(1,:),'-')==5)
            	formstr = 'yyyy-mm-dd HH:MM:SS';
          	case length(findstr(date(1,:),'-'))==2 & any(findstr(date(1,:),'-')==3) & str2double(date(1,1:2))<13
            	formstr = 'mm-dd-yyyy HH:MM:SS';
            otherwise
               	error('Date format is not recognized')
        end
        % Now convert the dates using the proper format
        jd = datenum(date,formstr);
    case date(1)>729756 & date(1)<737791
        % Julian day format
        jd = date;
    case date(1)>0 & date(1)<10000
        % Corie day format
        jd = date + 729023;
    otherwise
        error('Date format is not recognized')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varname = var_process(var)
% Determine which variable to extract
switch true
    case strcmp(var(1:2),'te')|strcmp(var(1:2),'Te')|strcmp(var(1:2),'TE')
        varname = 'temp.63';
    case strcmp(var(1:2),'sa')|strcmp(var(1:2),'Sa')|strcmp(var(1:2),'SA')
        varname = 'salt.63';
    case strcmp(var(1:2),'de')|strcmp(var(1:2),'De')|strcmp(var(1:2),'DE')
        varname = 'conc.63';
    case strcmp(var(1:2),'el')|strcmp(var(1:2),'El')|strcmp(var(1:2),'EL')
        varname = 'elev.61';
    case strcmp(var(1:3),'vel')|strcmp(var(1:3),'Vel')|strcmp(var(1:3),'VEL')
        varname = 'hvel.64';
    otherwise
        varname = var;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y] = pos_process(ypos,xpos)

% Ensure that lat's and lon's were entered correctly
if (abs(xpos(1))<90 & abs(ypos(1))<180)
	temp = xpos;
	xpos = ypos;
	ypos = temp;
 	clear temp
end
xpx = size(xpos,1);
xpy = size(xpos,2);
ypx = size(ypos,1);
ypy = size(ypos,2);
switch true
    case (xpx==1 & xpy==1) & (ypx>1 | ypy>1) % xpos only is a singleton
    	xpos = ones(size(ypos)).*xpos;
    case (ypx==1 & ypy==1) & (xpx>1 | xpy>1) % ypos only is a singleton
    	ypos = ones(size(xpos)).*ypos;
    case length(xpos)~=1 & length(ypos)~=1  % Build matrix of xy values 
    	[xpos,ypos] = meshgrid(xpos,ypos);
    case length(xpos)~=1 & length(xpos)==length(ypos)  % xy values are 1-to-1 
    	% No conversion is required
end
if ~any(abs(xpos)>180) & ~any(abs(ypos)>90)
    if xpos(1)>0
    	xpos = -xpos;
    end
    [x,y] = convll2m(xpos,ypos);
else
    x = xpos;
    y = ypos;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [db,dir] = db_process(db_name)

switch true
    case strcmp(db_name,'db26')|strcmp(db_name,'DB26')|strcmp(db_name,'26')
        db = '26';
        dir = '/home/workspace/ccalmr/hindcasts/';
    case strcmp(db_name,'db22')|strcmp(db_name,'DB22')|strcmp(db_name,'22')
        db = '22';
        dir = '/home/workspace/ccalmr/hindcasts/';
    case strcmp(db_name,'db16')|strcmp(db_name,'DB16')|strcmp(db_name,'16')
        db = '16';
        dir = '/home/workspace/ccalmr/hindcasts/';
    case strcmp(db_name,'db14')|strcmp(db_name,'DB14')|strcmp(db_name,'14')
        db = '14';
        dir = '/home/workspace/ccalmr/hindcasts/';
    case strcmp(db_name,'dev')|strcmp(db_name,'DEV')|strcmp(db_name,'Dev')
        db = 'dev';
        dir = '/home/workspace/ccalmr/hindcasts/';
    case strcmp(db_name,'f22')|strcmp(db_name,'F22')
        db = 'f22';
        dir = '/home/workspace/local0/forecasts/';
    case strcmp(db_name,'f26')|strcmp(db_name,'F26')
        db = 'f26';
        dir = '/home/workspace/local0/forecasts/';
    case any(findstr(db_name,'/')) | strcmp(db_name(1),'.')% This is a specific run
        dir = db_name;
        % Need to determine the start date of this run
        if exist([dir '/bctides.in'],'file')% Get info from the bctides.in file and pray it's correct!
            fid = fopen([dir '/bctides.in']);
            if fid<0
                error('Sorry, you need to ensure that a valid bctides.in file is in the requested run directory');
            end
            dstr = fgets(fid);
            db = datenum(dstr(1:end-4),'mm/dd/yyyy HH:MM:SS');
            fclose(fid);
        end
    otherwise
        error('Sorry, I don''t recognize that database source. Please check your input string.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z = depth_process(z)

% Check to make sure all depths are from the surface
if any(z<0) & any(z>0)
    error('All depths must be from the free surface (do not mix positive and negative values)')
end

if ~ischar(z(1))
    if any(z<0)
    	z = -z;
    end
    if size(z,2)>size(z,1)
        z = z';
    end
else
    switch true
        case strcmp(z(1:3),'Bot')|strcmp(z(1:3),'bot')
            z = -1;
        case strcmp(z(1:3),'Ave')|strcmp(z(1:3),'ave')|strcmp(z(1:3),'Mea')|strcmp(z(1:3),'mea')
            z = -2;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [jd,x,y,z] = redimensionalize(jd,x,y,z)

switch true
    case length(jd)==1&&length(x)==1&&length(z)~=1% CTD drop
        jd = ones(size(z)).*jd; x = ones(size(z)).*x;y = ones(size(z)).*y;
    case length(jd)==1&&length(x)~=1&&length(z)~=1% Vertical transect
        [x,fill] = meshgrid(x,z); [y,z] = meshgrid(y,z); jd = ones(size(x)).*jd;
    case length(jd)~=1&&length(x)==1&&length(z)~=1% Saturn1 style timeseries
        [jd,z] = meshgrid(jd,z); x = ones(size(jd)).*x; y = ones(size(jd)).*y;
    case length(jd)~=1&&length(x)~=1&&length(z)==1% Ship's TSG style plot
        z = ones(size(jd)).*z;
    case length(jd)~=1&&length(x)==1&&length(z)==1% Station timeseries
        x = ones(size(jd)).*x; y = ones(size(jd)).*y; y = ones(size(jd)).*y; z = ones(size(jd)).*z;
    case length(jd)==1&&length(z)==1&&size(x,1)>1|size(x,2)>1% Aerial view
        jd = ones(size(x)).*jd; z = ones(size(x)).*z;
    case length(jd)==length(z)&&length(z)==length(x)&&length(x)==length(y)% xyzt format
        jd = jd; x = x; y = y; z = z;
    otherwise
        error('The structure of the input arguments is not understood')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mode,t,x,y,z] = mode_process(jd,xpos,ypos,zpos)

switch true
    case zpos(1)==-1
        mode = 'xyt_bot';
        t = jd; x = xpos; y = ypos; z = zeros(size(jd));
    case zpos(1)==-2
        mode = 'xyt_ave';
        t = jd; x = xpos; y = ypos; z = zeros(size(jd)); 
    case zpos(1)==-3
        mode = 'vert_profile';
        t = jd; x = xpos; y = ypos; z = zeros(size(jd));
    case (length(xpos)==length(ypos)==length(zpos)==length(jd)) | ~ischar(zpos)
        mode = 'xyzt';
        t = jd; x = xpos; y = ypos; z = zpos;        
    otherwise
        error('The structure of the input arguments is not understood')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fid = build_input(t,x,y,z,var,db,dir,mode)
switch true
    case strcmp(mode,'xyzt') | strcmp(mode,'xyt_ave') | strcmp(mode,'xyt_bot');
        s = rem(t,1).*86400;s(s==0)=86400;
        n = length(t);
        % Determine which file to access
        if ischar(db)% This is a standard database
            dv = datevec(t(1)); year=dv(1); yr = num2str(year); sol_day = floor(t(1)-datenum(year,1,1)+1);
            week = floor((sol_day-1)/7)+1; wk = num2str(week,'%.2i'); corie_start = floor(t(1))-729023; 
            iter = mod(sol_day-1,7)+1;
            file = [num2str(iter) '_' var];
            fildir = [dir yr '-' wk '-' db '/run/' file];
            if ~exist(fildir,'file')% Then it may be a forecast directory
                % Calculate solar day
                file = ['2_' var];
                fildir = [dir db '/' yr '-' num2str(sol_day,'%.3i') '/run/' file];
            end
        else
            jd0 = datenum(db);% In this case, 'db' is the model start time in julian days
            jd1 = datenum(t(1));
            corie_start = floor(t(1))-729023; 
            iter = floor((jd1-jd0)+1);
            file = [num2str(iter) '_' var];
            fildir = [dir '/outputs/' file];
        end
        if ~exist(fildir,'file')
            error('Model results unavailable for that time/database combination')
        end
        block = [(1:n)' x y z s];
        fid = fopen('modext_input.dat','w');
        fprintf(fid,'%s\n',fildir);
        fprintf(fid,'%i\n',corie_start);
        fprintf(fid,'%i\n',n);
        fprintf(fid,'%i %f %f %f %f\n',block');
        fclose(fid);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y]=convll2m(lon,lat) 

outsize = size(lon);
lon = lon(:);
lat = lat(:);
data = [(1:length(lon));lon';lat';ones(1,length(lat))];
fid = fopen('temp.dat','w');
fprintf(fid,'%s\n','***');
fprintf(fid,'%i\n',length(lon));
fprintf(fid,'%12i %12.3f %12.3f %i\n',data);
fclose(fid);

!LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH  /home/users/cseaton/bin/spcs2ll_bp  -input   temp.dat    -output    temp2.dat    -ll2spcs  > /dev/null 2>&1
fid = fopen('temp2.dat','r');
fgets(fid);
fgets(fid);
out = fscanf(fid,'%i %f %f %f\n',[4 inf]);
x = out(2,:);
y = out(3,:);
xt = zeros(outsize);
yt = xt;
xt(:) = x(:);
yt(:) = y(:);
x = xt;
y = yt;
!rm -f temp.dat temp2.dat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = fullfile_extract(db_name,date,zpos)
% This function needs additional functionality. Code to extract values at 
% constant depths needs to be developed. 

out = melioext(db_name,date);
if ~isstr(zpos)
    out = out; % Interpolate data to a particular depth
else
    switch true
        case (strcmp(zpos,'All')|strcmp(zpos,'all')) | str2num(db_name(end))<3 % Take everything returned
            out = out;
        case (strcmp(zpos(1:3),'Bot')|strcmp(zpos(1:3),'bot')) & str2num(db_name(end))>2 % Take bottom values
            for i = 1:length(out.data)
                if isfield(out.data,'val')
                    out.data(i).val = out.data(i).val(:,end);
                else
                    out.data(i).u = out.data(i).u(:,end);
                    out.data(i).v = out.data(i).v(:,end);
                end
                out.data(i).z = out.data(i).z(:,end);
            end
        case (strcmp(zpos(1:3),'Sur')|strcmp(zpos(1:3),'sur')) & str2num(db_name(end))>2 % Take surface values
            for i = 1:length(out.data)
                if isfield(out.data,'val')
                    out.data(i).val = out.data(i).val(:,1);
                else
                    out.data(i).u = out.data(i).u(:,1);
                    out.data(i).v = out.data(i).v(:,1);
                end
                out.data(i).z = out.data(i).z(:,1);
            end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
