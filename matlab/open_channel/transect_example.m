% Example use of transect routines
%
path(path,'/usr/local/cmop/matlab/cmop/m-elio');
%path(path,'../');
tr.hgrid=gr_readHGrid('n_channel.bp');
gr.hgrid=gr_readHGrid('hgrid.gr3');
h=sz_readHeader('outputs/4_salt.63');
he=sz_readHeader('outputs/4_elev.61');
gr.vgrid=h.vgrid;

% compute transect
[ob]= ob_ini_fromTrasect(gr, 'n_channel.bp');
trLen = cumsum(sqrt((ob.xy.x(2:end)-ob.xy.x(1:end-1)).^2+(ob.xy.y(2:end)-ob.xy.y(1:end-1)).^2));
trLen = [0; trLen];

% Read timestep of data
ds=sz_readTimeStep(h,1);
dsd=map_sz2hts(h, ds);

s=ob.xy.H*double(dsd);

% Construct vertical grid
% Read timestep of elevations
de=sz_readTimeStep(he,1);

e=ob.xy.H*double(de);

dp=ob.xy.H*gr.hgrid.depth;

sz=sz_computeZlevels(dp,e,gr.vgrid);

h=pcolor(trLen, sz', s');
set(h,'EdgeColor','none','FaceColor','interp');
axis([-inf inf -30 2])
caxis([0 32]); hold on;
