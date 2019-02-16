function dir=safeDir(dirIn,dirFlag)
%
% dir=safeDir(dirIn)
% function to replace file seperator with the correct os dependent
% seperator (i.e. "/" windows, "\" linux/unix...)
%
% if optional dirFlag=1 then function will chec for ending filesep and add
% if needed.
%
% benk 2011

if nargin<2, dirFlag=0; end

osSep=filesep;
switch osSep
    case '/'
        otherSep='\';
    case '\'
        otherSep='/';
end

dir=strrep(dirIn,otherSep,osSep);
if dirFlag
    if dir(end)~=filesep
        dir=[dir filesep];
    end
end


