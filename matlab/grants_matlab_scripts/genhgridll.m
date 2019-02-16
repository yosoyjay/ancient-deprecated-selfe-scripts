fg = hgrid2fg('hgrid.gr3');
fid = fopen('hgrid.ll','w');
fprintf(fid,'%s\n','Hgrid in lat/lon');
fprintf(fid,'%i %i\n',[length(fg.e) length(fg.x)]);
[lon,lat] = convm2ll(fg.x,fg.y);
block = [1:length(fg.x);lon';lat';fg.z'];
fprintf(fid,'%i %f %f %f\n',block);
block = [1:length(fg.e);ones(1,length(fg.e)).*3;fg.e'];
fprintf(fid,'%i %i %i %i %i\n',block);
fclose(fid);
status = exist('hgrid.ll','file');
if status~=2
    error('hgrid.ll generation unsuccessful')
end