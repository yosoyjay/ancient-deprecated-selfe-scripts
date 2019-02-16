function [wind,wfac] = gengr3(hgrid_dir,new_hgrid)

fg1 = hgrid2fg([hgrid_dir '/hgrid.gr3']);
fg2 = hgrid2fg(new_hgrid);
x1=fg1.x; y1=fg1.y; x2=fg2.x; y2=fg2.y; 
nodes = length(x2);
elems = length(fg2.e);
wind = ones(nodes,3);
wfac = wind;
h = waitbar(0,'Creating interpolation weighting data...');
for i = 1:nodes
    dist = sqrt((x1-x2(i)).^2 + (y1-y2(i)).^2);
    [dist,ind] = sort(dist);
    if dist(1)>1
        wind(i,:) = ind(1:3);
        inv_d = 1./dist(1:3);
        wfac(i,:) = inv_d ./ sum(inv_d);
    else
        wind(i,:) = ind(1:3);
        wfac(i,:) = [1 0 0];
    end
    waitbar(i/nodes,h)
end
close(h);
save temp
n = [(1:nodes)' x2 y2 fg2.z];
e = [(1:elems)' ones(elems,1).*3 fg2.e];

filename = {'diffmax.gr3' 'diffmin.gr3' 'drag.gr3' 'estuary.gr3' 'interpol.gr3' 's_nudge.gr3' 'watertype.gr3' 'albedo.gr3' 'windrot_geo2proj.gr3' 'xlsc.gr3' 't_nudge.gr3'};
for i = 1:length(filename)
	fid3 = fopen([hgrid_dir '/' filename{i}],'r');
    fgets(fid3);
    out = fscanf(fid3,'%i %i\n',2);
    data = fscanf(fid3,'%i %f %f %f\n',[4,out(2)]);
    data = data(4,:)';
    if sum(mod(data,1))<.0001
        fstr = '%16i %16.4f %16.4f %16i\n';
        new = round(sum(data(wind') .* wfac')');
    else
        fstr = '%16i %16.4f %16.4f %16.4f\n';
        new = sum(data(wind') .* wfac')';
    end
    n(:,4) = new;
	if ~exist(filename{i},'file')
        fid = fopen(filename{i},'w');
    else
    	error('Please discard current GR3 files');
    end
	fprintf(fid,'%s\n',filename{i});
 	fprintf(fid,'%16i %16i\n',[elems nodes]);
  	fprintf(fid,fstr,n');
  	fprintf(fid,'%16i %16i %16i %16i %16i\n',e');
end

[lon,lat] = convm2ll(x2,y2);
fid = fopen('hgrid.ll','w');
fprintf(fid,'%s\n','hgrid.ll');
fprintf(fid,'%i %i\n',[elems nodes]);
fprintf(fid,'%i %f %f %f\n',[(1:nodes);lon';lat';fg2.z']);
fprintf(fid,'%i %i %i %i %i\n',[(1:elems);ones(1,elems).*3;fg2.e']);
fclose(fid);
!rm temp.mat
