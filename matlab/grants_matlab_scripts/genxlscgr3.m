function status = genxlscgr3(val,method)
% 

fg = hgrid2fg('hgrid.gr3');
fid = fopen('xlsc.gr3','w');
fprintf(fid,'%s\n','xlsc values');
fprintf(fid,'%i %i\n',[length(fg.e) length(fg.x)]);

switch method
    case 'constant'
        block = [1:length(fg.x);fg.x';fg.y';ones(1,length(fg.x)).*val];
        fprintf(fid,'%i %f %f %f\n',block);
    otherwise
        error('Didn''t recognize that method...');
end

block = [1:length(fg.e);ones(1,length(fg.e)).*3;fg.e'];
fprintf(fid,'%i %i %i %i %i\n',block);
fclose(fid);
status = exist('watertype.gr3','file');