close all
path(path,'./m-elio');
%--Coastline Information
load 'c17843.dat';
c1=c17843;
x1=c1(:,1);y1=c1(:,2);

% prepare fields - salt & veolocity
gr.hgrid=gr_readHGrid('hgrid.gr3');
gr2.hgrid=gr_readHGrid('hgrid.ll');
fg = hgrid2fg('hgrid.ll');

h=sz_readHeader('3_salt.63');
gr.vgrid=h.vgrid;
hv=sz_readHeader('3_hvel.64');

% set some default values
aa = [47.472609 63.592969 19.899843 30.46383]; % map region

it = 1
vbot = 18
[d ts]=sz_readTimeStep(h,it); %step 64 is 1600 hrs PST (= 0000 GMT)
[u] = map_sz2hts(h,d,1);
[s0] = u(:,end); % End is the surface layer
[s1] = u(:,vbot); % vbot is the bottom layer
s0(s0<0) = NaN;
s1(s1<0) = NaN;
sal=s0;
salb=s1;

dlat0 = 10; % spacing for regridded velocity (degrees of latitude)
[LON0,LAT0] = meshgrid([aa(1):dlat0:aa(2)],[aa(3):dlat0:aa(4)]);

sg0 = griddata(fg.x,fg.y,s0,LON0,LAT0,'cubic');
sg1 = griddata(fg.x,fg.y,s1,LON0,LAT0,'cubic');

[d ts]=sz_readTimeStep(hv,it); %step 64 is 1600 hrs PST (= 0000 GMT)
u = map_sz2hts(hv,d(:,1),1);
v = map_sz2hts(hv,d(:,2),1);
u = u(:,end);
v = v(:,end);

utop = griddata(fg.x,fg.y,u,LON0,LAT0,'cubic');
vtop = griddata(fg.x,fg.y,v,LON0,LAT0,'cubic');
utop0=utop;vtop0=vtop;

%dar = [1/cos(pi*46/180) 1 1];

z0 = griddata(fg.x,fg.y,-gr.hgrid.depth,LON0,LAT0,'cubic');

dlat2 = 100;
[LON2,LAT2] = meshgrid([aa(1):dlat2:aa(2)],[aa(3):dlat2:aa(4)]);
dlat = 0.1; % spacing for regridded velocity (degrees of latitude)
[LON,LAT] = meshgrid([aa(1):dlat:aa(2)],[aa(3):dlat:aa(4)]);


uu = interp2(LON0,LAT0,utop0,LON,LAT);
vv = interp2(LON0,LAT0,vtop0,LON,LAT);
uu2 = interp2(LON0,LAT0,utop,LON2,LAT2);
vv2 = interp2(LON0,LAT0,vtop,LON2,LAT2);

slo = 35; shi = 40;
caxis([slo shi]);
colorbar('north')
title('(a) Surface Salinity');
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');

gr_plot2(gr2.hgrid,sal,[slo shi]);

% add velocity vectors
ufact = 1.0;
quiver(LON,LAT,ufact*uu,ufact*vv,0,'k');

