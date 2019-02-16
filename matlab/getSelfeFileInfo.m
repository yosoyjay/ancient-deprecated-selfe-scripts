function dat=getSelfeFileInfo(dirIn)

olddir=cd;
cd(safeDir(dirIn))
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
    dat.filetype=[dat.filetype{i}]'; % file types
    dat.filenums=str2double([unique([dat.filenums{:}])']);  %time index
    
    dat.nVariables=length(dat.varNames);    
    dat.nFiles=length(dat.filenums);    
else
    dat=[];
end
    
    
cd(olddir);