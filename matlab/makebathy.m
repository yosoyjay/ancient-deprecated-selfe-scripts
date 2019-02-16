fg = hgrid2fg('shared_data/200ea_hgrid.gr3');
fgsub = hgrid2fg('shared_data/subsidence_hgrid.gr3');
fgpdev = hgrid2fg('shared_data/predev_hgrid.gr3');

narr = fg.z;
narr = narr-max(narr);
narr(narr<-50) = -50;
narr = narr.*-64/50;
narr = ceil(narr);
narr(narr==0) = 1;
narr = 65-narr;
jarr = jet;
c = jarr(narr,:);
set(gcf,'Position',[179         221        1156         603])
set(gca,'Position',[.01 .02 .98 .96])
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])
set(gca,'Color',[.15 .15 .15])
cfloor = patch('Vertices',[fg.x fg.y fg.z],'Faces',fg.e,'FaceVertexCData' ...
,c,'FaceColor','interp','EdgeColor','none', ...
'SpecularStrength',0,'DiffuseStrength',.7);
L = light;
lighting phong
lightangle(L,180,60)
set(gca,'XLim',[328980 403570])
set(gca,'YLim',[265150 308850])
set(gca,'ZLim',[-100 100])
material dull

figure
narr = fgsub.z;
narr = narr-max(narr);
narr(narr<-50) = -50;
narr = narr.*-64/50;
narr = ceil(narr);
narr(narr==0) = 1;
narr = 65-narr;
c = jarr(narr,:);
set(gcf,'Position',[179         221        1156         603])
set(gca,'Position',[.01 .02 .98 .96])
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])
set(gca,'Color',[.15 .15 .15])
cfloor = patch('Vertices',[fgsub.x fgsub.y fgsub.z],'Faces',fgsub.e,'FaceVertexCData' ...
,c,'FaceColor','interp','EdgeColor','none', ...
'SpecularStrength',0,'DiffuseStrength',.7);
L = light;
lighting phong
lightangle(L,180,60)
set(gca,'XLim',[328980 403570])
set(gca,'YLim',[265150 308850])
set(gca,'ZLim',[-100 100])
material dull

figure
narr = fgpdev.z;
narr = narr-max(narr);
narr(narr<-50) = -50;
narr = narr.*-64/50;
narr = ceil(narr);
narr(narr==0) = 1;
narr = 65-narr;
c = jarr(narr,:);
set(gcf,'Position',[179         221        1156         603])
set(gca,'Position',[.01 .02 .98 .96])
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])
set(gca,'Color',[.15 .15 .15])
cfloor = patch('Vertices',[fgpdev.x fgpdev.y fgpdev.z],'Faces',fgpdev.e,'FaceVertexCData' ...
,c,'FaceColor','interp','EdgeColor','none', ...
'SpecularStrength',0,'DiffuseStrength',.7);
L = light;
lighting phong
lightangle(L,180,60)
set(gca,'XLim',[328980 403570])
set(gca,'YLim',[265150 308850])
set(gca,'ZLim',[-100 100])
material dull