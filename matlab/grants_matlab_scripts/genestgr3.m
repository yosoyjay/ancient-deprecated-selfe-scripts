s = [-124.07785735 46.2521884]; % Ocean side of initial estuarine salinity gradient
f = [-123.75586705 46.21700909];% River side of initial estuarine salinity gradient
nj = [-124.08394402 46.26533149];% Point position of North Jetty
sj = [-124.07264769 46.23380919];% Point position of South Jetty
fg = hgrid2fg('hgrid.gr3');
fg.bndpth(isnan(fg.bndpth)) = [];

% Find polygon that encircles the estuary
[fg.lon,fg.lat] = convm2ll(fg.x,fg.y);
% North Jetty
diff = sqrt((fg.lon(fg.bndpth)-nj(1)).^2 + (fg.lat(fg.bndpth)-nj(2)).^2);
njnode = find(diff == min(diff)); njnode = njnode(1);
% South Jetty
diff = sqrt((fg.lon(fg.bndpth)-sj(1)).^2 + (fg.lat(fg.bndpth)-sj(2)).^2);
sjnode = find(diff == min(diff)); sjnode = sjnode(1);
if sjnode<njnode
    estpoly = [fg.lon(fg.bndpth(sjnode:njnode)) fg.lat(fg.bndpth(sjnode:njnode))];
else
    estpoly = [fg.lon(fg.bndpth(njnode:sjnode)) fg.lat(fg.bndpth(njnode:sjnode))];
end
estpoly(end+1,:) = estpoly(1,:);
est = inpolygon(fg.lon,fg.lat,estpoly(:,1),estpoly(:,2));
% Set s and f indices
diff = sqrt((fg.lon-s(1)).^2 + (fg.lat-s(2)).^2);
sn = find(diff==min(diff)); sn = sn(1);
diff = sqrt((fg.lon-f(1)).^2 + (fg.lat-f(2)).^2);
fn = find(diff==min(diff)); fn = fn(1);
est = double(est);
est(sn) = -1;
est(fn) = -2;

fid = fopen('estuary.gr3','w');
fprintf(fid,'%s\n','Estuary definitions for initialization');
fprintf(fid,'%i %i\n',[length(fg.e) length(fg.x)]);
block = [1:length(fg.x);fg.x';fg.y';est'];
fprintf(fid,'%16i %16.4f %16.4f %16i\n',block);
block = [1:length(fg.e);ones(1,length(fg.e)).*3;fg.e'];
fprintf(fid,'%16i %16i %16i %16i %16i\n',block);
fclose(fid);
status = exist('estuary.gr3','file');
if status~=2
    error('estuary.gr3 generation unsuccessful')
end
