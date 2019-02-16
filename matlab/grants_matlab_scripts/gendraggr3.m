function status = gendraggr3(drag,method)
% 

fg = hgrid2fg('hgrid.gr3');
fid = fopen('drag.gr3','w');
fprintf(fid,'%s\n',['Drag values using ' method ' method']);
fprintf(fid,'%i %i\n',[length(fg.e) length(fg.x)]);

switch method
    case 'constant'
        if length(drag)>1
            error('For the "constant" method, drag must be scalar')
        end
        block = [1:length(fg.x);fg.x';fg.y';ones(1,length(fg.x)).*drag];
        fprintf(fid,'%i %f %f %f\n',block);
    case 'depth' % This methods attempts to reproduce whatever Joe Cho was doing in his runs???
        z = fg.z;
        z(fg.z>0) = 0;
        mind = drag(1); maxd = drag(2);
        drag = exp((z-45)./8) +mind;
        drag(drag>maxd) = maxd;
        block = [1:length(fg.x);fg.x';fg.y';drag'];
        fprintf(fid,'%i %f %f %f\n',block);
    otherwise
        error('Didn''t recognize that method...');
end

block = [1:length(fg.e);ones(1,length(fg.e)).*3;fg.e'];
fprintf(fid,'%i %i %i %i %i\n',block);
fclose(fid);
status = exist('drag.gr3','file');