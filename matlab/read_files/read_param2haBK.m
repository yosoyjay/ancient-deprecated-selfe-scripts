function param=read_param2haBK(file)
%
%function to read in freeformat param files
%
%ben knight 2011

fid=fopen(file,'r');
count=0
while 1
    count=count+1;
    s=fgetl(fid);
    if ~ischar(s), break, end
    s=strtrim(s);
    if ~isempty(s)
    if s==-1, break, end;
    if ~strcmp(s(1),'!')
        s=strrep(s,'!','; %');
        if isempty(regexp(s,'\w[.]\w*\s*[=]'))   
            try 
                eval(['param.' s ';']);
            catch
                warning(['Unreconized data in param file, line: ' num2str(count)])
            end
        end
    end
    end
    end
end

