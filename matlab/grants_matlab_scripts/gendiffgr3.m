function status = gendiffgr3(min,max,method)
% 

fg = hgrid2fg('hgrid.gr3');
fid1 = fopen('diffmin.gr3','w');
fid2 = fopen('diffmax.gr3','w');
fprintf(fid1,'%s\n','Minimum diffusion values');
fprintf(fid2,'%s\n','Maximum diffusion values');
fprintf(fid1,'%i %i\n',[length(fg.e) length(fg.x)]);
fprintf(fid2,'%i %i\n',[length(fg.e) length(fg.x)]);

switch method
    case 'constant'
        block1 = [1:length(fg.x);fg.x';fg.y';ones(1,length(fg.x)).*min];
        block2 = [1:length(fg.x);fg.x';fg.y';ones(1,length(fg.x)).*max];
        fprintf(fid1,'%i %f %f %f\n',block1);
        fprintf(fid2,'%i %f %f %f\n',block2);
    otherwise
        error('Didn''t recognize that method...');
end

block = [1:length(fg.e);ones(1,length(fg.e)).*3;fg.e'];
fprintf(fid1,'%i %i %i %i %i\n',block);
fprintf(fid2,'%i %i %i %i %i\n',block);
fclose(fid1);
fclose(fid2);
status(1) = exist('diffmin.gr3','file');
status(2) = exist('diffmax.gr3','file');