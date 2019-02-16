%%%%%%% generic script to produce images for on-line ptrack server %%%%%%%

%%%% hard-coded args: %%%%
fparams = 'image_params.dat';

%%% params is a struct including fields: trackid, itype, fname
params = read_ptrack_images_params(fparams);

pathdataxy = partpath('particle.pth');

pathdata = partpath_ll('particle.ll');

obs=[];
obs_time=[];
clustsz=[];


plottype = 'Simulation';
% if (strcmp(params.itype(1:4), 'base') || strcmp(params.itype(1:4), 'basi')) 
%     plottype = 'Simulation';
if (length(params.itype)>=7 & strcmp(params.itype(1:7), 'drifter'))
    [tmp, lon, lat, dnum] = read_drifter('drifter.dat');
    obs=[lon', lat'];
    obs_time = dnum;
end

include_pt = [-124.01388, 46.248974];   %include this point to make sure there is geographical context

% hpath = plotpath(plottype,relshape,figure_title,clustsz);
hpath = plotpath_ll(plottype, size(pathdata.x,1), pathdata.time(1,1),...
                   [params.trackid ' particle track'], clustsz,include_pt,obs,obs_time);
% sdat = hgexport('factorystyle');
% hgexport(gcf,'ptrack_plot.eps',sdat,'format','eps');
set(gcf, 'PaperPositionMode', 'auto'); 
print('-dpng', [params.fname '.png']);

