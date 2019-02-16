fg = hgrid2fg('hgrid.gr3');
fid = fopen('windrot_geo2proj.gr3', 'w');
fprintf(fid,'%s\n','Rotation correction factors for wind data provided in lat/lon format');
fprintf(fid,'%i %i\n',[length(fg.e) length(fg.x)]);
int = ones(size(fg.z));
int(fg.z>sz_trans) = 2;
block = [1:length(fg.x);fg.x';fg.y';int'];
fprintf(fid,'%i %f %f %i\n',block);
block = [1:length(fg.e);ones(1,length(fg.e)).*3;fg.e'];
fprintf(fid,'%i %i %i %i %i\n',block);
fclose(fid);
status = exist('interpol.gr3','file');
if status~=2
    error('interpol.gr3 generation unsuccessful')
end
