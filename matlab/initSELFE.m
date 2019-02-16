function theResult = initSELFE(DirIn)
% function to initialise selfe object with 

if nargin < 1, help(mfilename), return, end
DirIn=safeDir(DirIn);
DirIn=strtrim(DirIn);
osSep=filesep;
if ~strcmp(DirIn(end),osSep)
    DirIn=[DirIn osSep];
end

status=exist(DirIn,'dir');
if status >= 0
    olddir=cd;
    cd(DirIn)
    
    %static data
    result.param_file='param.in';
    result.grid_file='hgrid.gr3';
    result.vgrid_file='vgrid.in';
    result.dateFormat='dd/mm/yyyy HH:MM'
    
    %get dir information
   result=getSelfeFileInfo(DirIn,result);
   
   %read in key grid/time information
   result=getGridInfo(DirIn,result);
   
   cd(olddir);
end
if nargout > 0, theResult = result; end


%% helper functions
function dat=getSelfeFileInfo(dirIn,dat)
%get info on output files

varnames=dir('*.6*');
if ~isempty(varnames)
    dat.SELFE_Dir=safeDir(dirIn);
    out={varnames.name}';
    dat.filetype=regexp(out,'[.]+\d+','match');
    dat.filenums=regexp(out,'^\d+','match');
    dat.varNames=regexp(out,'[_]\w+[_]?\d*[.]','match');

    %get output file details
    [dat.varNames i j]=unique([dat.varNames{:}]); %variable names
    dat.varNames=dat.varNames';
    for k=1:length(dat.varNames)
        dat.varNames{k}=dat.varNames{k}(2:(end-1));
    end
    dat.filetype=[dat.filetype{i}]'; % file types
    dat.filenums=sort(str2double([unique([dat.filenums{:}])']));  %time index
    
    dat.nVariables=length(dat.varNames);    
    dat.nFiles=length(dat.filenums);    
else
    dat=[];
end



function dat=getGridInfo(dirIn,dat)
% get grid and time info from files in directory

%open files
dat.param=read_param2haBK(dat.param_file);
dat.hgrid=gr_readHGrid(dat.grid_file);
dat.vgrid=read_vgridbk(dat.vgrid_file);
dat.vgrid.sLevels=dat.vgrid.Scoord;
dat.vgrid.h0=dat.param.h0;

%get start date from start of param file
fid=fopen(dat.param_file);
strline=fgetl(fid);
while ~isdate(strline)
    strline=fgetl(fid);
end
fclose(fid);
dat.param.startdate=datenum(strline(2:end),dat.dateFormat); %note use (2:end) to skip leading "!"

%get time information
% Time steps between file output
dat.nspool=dat.param.nspool; %60;
%The number iteration per file 
dat.iter=dat.param.ihfskip/dat.param.nspool;

% create a time index
dat.tstartFiles=dat.iter.*(dat.filenums-1).*dat.nspool+1;  % in timestep units

dat.time=zeros(length(dat.filenums)*dat.iter,3); %time,file,idx1,idx2
for i=1:length(dat.tstartFiles)
    idx=(((i-1)*dat.iter)+1):((i*dat.iter));
    dat.time(idx,1)= dat.tstartFiles(i)+(1:dat.iter)*dat.nspool;  %time in timesteps of outputs, file index
    dat.time(idx,2)= dat.filenums(i); %filenum
    dat.time(idx,3)=1:dat.iter; %idx1
end
dat.time(:,4)=1:(length(dat.filenums)*dat.iter); %idx2
dat.time(:,1)=dat.time(:,1)*dat.param.dt;  %convert to seconds

dat.datetime= dat.param.startdate+dat.time(:,1)/24/60/60; %matlab date time
dat.maxi_iter=max(dat.time);




