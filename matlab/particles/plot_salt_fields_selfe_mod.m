% plot_salt_fields_selfe.m  11/19/2008 Parker MacCready
%
% This plots a figure of model surface salinity, a salinity section, and an
% observed salinity section
%
% The intention is to hand this to Antonio's group and they will add in
% panels from their model (SELFE)
%

%clear
close all

path(path,'./m-elio');

%--Coastline Information
load '/home/workspace/ccalmr45/choj/for_grant/barb_sal_maps/c28142.dat';
c1=c28142;
x1=c1(:,1);y1=c1(:,2);

%--Loading CR-line stations
%load '/home/workspace/ccalmr45/choj/for_grant/barb_sal_maps/crline_aug31.in';
%c2=crline_aug31;
%--Loading GH-line stations
%load '/home/workspace/ccalmr45/choj/for_grant/barb_sal_maps/ghline_sep1.in';
%c2=ghline_sep1;
%--Loading LP-line stations

load '/home/workspace/ccalmr45/choj/for_grant/barb_sal_maps/lpline_sep2.in';
c2=lpline_sep2;
x2=-c2(:,2);y2=c2(:,1);

load '/home/workspace/ccalmr45/choj/for_grant/barb_sal_maps/mission_28.pts';
load '/home/workspace/ccalmr45/choj/for_grant/barb_sal_maps/mission_17.pts';
c3=mission_28;
x3=c3(:,1);y3=c3(:,2);

c4=mission_17;
x4=c4(:,1);y4=c4(:,2);

base_out_dir ='/home/workspace/ccalmr45/choj/for_grant/barb_sal_maps/';
start_day = datenum(2010, 1, 1, 0, 0, 0);

%-- This section will be changed for each day and each forecast.
%iopt=2;
curr_run_links_dir = '/home/workspace/local0/forecasts';
if iopt==1
  modnames={'f22'};
end
if iopt==2
  modnames={'dev'};
end
if iopt==3
  modnames={'f14'};
end
if iopt==4
  modnames={'db16'};
end

%yr=2009;imn=11;idy=20; % variables: today's date
% ist=start; intv=interval; ied=end; time loop: it=1:1:96--> every 15min; it=4:4:96--> hourly
% ik=1 for tomorrow; 2 for the next day
ist=1;intv=1;ied=96; ik=1;  

%Today
if ik==1
  fn = sprintf('%s-sal-%d-%02d-%02d.avi', modnames{1}, yr, imn, idy); 
  zz1=1;zz2=1;zz3=1;fnm='2_'; % for tomorrow
  aviobj=avifile(fn,'fps',0.5);
  %aviobj=avifile('ani4antonio_proposal.avi','fps',1);
end

% Next day...
% if ik==2
%    aviobj=avifile('ani_sal4oct10_fca200nb_estuary.avi','fps',1);
%    zz1=2;zz2=0;zz3=2;fnm='3_'; % for the day after tomorrow
% end

% Source of data (f22, dev, f14, db16).... vbot - vertical bottom
if iopt<=3
  vbot=18; 	% vbot=18 for fca200nb and dev;
else
  vbot=1; 	% for fdb21 and db16
end
%-------------------------------------------------------------

curr_day=datenum(yr, imn, idy, 0, 0, 0);
doy=curr_day-start_day;

for imod = 1:length(modnames)
model_day = sprintf('%d-%03d', yr, doy+zz2); 
model_day = sprintf('%d-%03d', yr, doy); 
disp(model_day);
model_day = 'today'; 
disp(model_day);
eval(['!ln -sf ' curr_run_links_dir '/' modnames{imod} '/' model_day '/run/hgrid.gr3 hgrid.gr3'])
eval(['!ln -sf ' curr_run_links_dir '/' modnames{imod} '/' model_day '/run/hgrid.ll hgrid.ll'])
eval(['!ln -sf ' curr_run_links_dir '/' modnames{imod} '/' model_day '/run/2_hvel.64 2_hvel.64'])
eval(['!ln -sf ' curr_run_links_dir '/' modnames{imod} '/' model_day '/run/2_salt.63 2_salt.63'])
eval(['!ln -sf ' curr_run_links_dir '/' modnames{imod} '/' model_day '/run/2_elev.61 2_elev.61'])
eval(['!ln -sf ' curr_run_links_dir '/' modnames{imod} '/' model_day '/run/3_hvel.64 3_hvel.64'])
eval(['!ln -sf ' curr_run_links_dir '/' modnames{imod} '/' model_day '/run/3_salt.63 3_salt.63'])
eval(['!ln -sf ' curr_run_links_dir '/' modnames{imod} '/' model_day '/run/3_elev.61 3_elev.61'])

s_list = {fnm,fnm,fnm}; 		% 2_ for Forecast Day 2; 3_ for forecast day 3
t_list = {'ghline2.bp','crline.bp','cmline.bp'};
sect_list = {'GH','CR','CM'};
%slat =[46.5,46.167,45.5];
%slonw=[-124.5,-125,-125];
%slone=[-124.05,-124.077,-123.95];
slat =[47,46.167,45.5]; 		%--47.91667,
slonw=[-125.4,-125,-125]; 		%---124.742,
slone=[-124.24,-124.077,-123.95]; 	%---125.928,
scalea={[0.6 0.7 0.37 0.2],[0.6 0.4 0.37 0.2],[0.6 0.1 0.37 0.2]};
nnn=0;
idy=idy+zz3; 

% prepare fields - salt & veolocity
gr.hgrid=gr_readHGrid('hgrid.gr3');
gr2.hgrid=gr_readHGrid('hgrid.ll');
fg = hgrid2fg('hgrid.ll');
h=sz_readHeader([s_list{3} 'salt.63']);
gr.vgrid=h.vgrid;
hv=sz_readHeader([s_list{3} 'hvel.64']);

% set some default values for the map
%   aa = [-125 -123.6 45 47]; 		% map region
    aa = [-125.5 -123.4 45 48.2]; 	% map region
    bb = [-124.5 -124 45.8 46.6]; 	% map zoom-in region
    ee = [-124.1 -123.3 46.1 46.4]; 	% map zoom-in region
    dar = [1/cos(pi*46/180) 1 1]; 	% Cartesian scaling
    clevs = [-200 -200]; 		% bathymetry contours
%   slo = 10; shi = 34; 		% salinity color range
    slo = 10; shi = 34; 		% salinity color range
    sclevs = [30:.5:34]; 		% salinity contour levels
    sclevs2 = [10:2:20]; 		% salinity contour levels
    fs1 = 14; 				% fontsize
    dlat0 = 0.001; 			% spacing for regridded velocity (degrees of latitude)

    [LON0,LAT0] = meshgrid([aa(1):dlat0*dar(1):aa(2)],[aa(3):dlat0:aa(4)]);
    z0 = griddata(fg.x,fg.y,-gr.hgrid.depth,LON0,LAT0,'cubic');
%   dlat2 = 0.05; 			% spacing for regridded velocity (degrees of latitude)
    dlat2 = 0.025;
    [LON2,LAT2] = meshgrid([aa(1):dlat2*dar(1):aa(2)],[aa(3):dlat2:aa(4)]);
    dlat = 0.1; 			% spacing for regridded velocity (degrees of latitude)
    [LON,LAT] = meshgrid([aa(1):dlat*dar(1):aa(2)],[aa(3):dlat:aa(4)]);

for it=ist:intv:ied  			% time loop: it=1:96--> every 15min; it=4:4:96--> hourly 
    nnn=nnn+1
    figure(100);
    if iopt<=3
       set(gcf,'position',[10 100 1400 750],'Color','w');
    else 
       set(gcf,'Position',[10 150 700 750],'Color','w');
    end
    ctime = [datestr(curr_day + (it * 900.0 / 86400.0), 'yyyy-mm-dd HH:MM:SS') ' PST'];

    set(0,'defaultaxesfontsize',fs1);
    set(0,'defaulttextfontsize',fs1);
    set(0,'defaultaxesfontname','Times');
    set(0,'defaulttextfontname','Times');
    colormap(jet(256)); % make a steppy colormap

    % SELFE Surface salinity map with velocity vectors
    [d ts]=sz_readTimeStep(h,it); 	%step 64 is 1600 hrs PST (= 0000 GMT)
    [u] = map_sz2hts(h,d,1);
    [s0] = u(:,end); 			% End is the surface layer
    [s1] = u(:,vbot); 			% vbot is the bottom layer
    s0(s0<0) = NaN;
    s1(s1<0) = NaN;
    sal=s0;
    salb=s1;
	
	%  For some data source
    if iopt<=2
	    s0 = griddata(fg.x,fg.y,s0,LON0,LAT0,'cubic');
	    s1 = griddata(fg.x,fg.y,s1,LON0,LAT0,'cubic');
	    s0(z0>-5) = NaN;
	    s0(LON0 > -124.25 & LAT0 > 47) = NaN;
    end

    % prepare fields - velocity
    [d ts]=sz_readTimeStep(hv,it); 	%step 64 is 1600 hrs PST (= 0000 GMT)
    u = map_sz2hts(hv,d(:,1),1);
    v = map_sz2hts(hv,d(:,2),1);
    u = u(:,end);
    v = v(:,end);
    utop = griddata(fg.x,fg.y,u,LON0,LAT0,'cubic');
    vtop = griddata(fg.x,fg.y,v,LON0,LAT0,'cubic');
    utop0=utop;vtop0=vtop;
    utop0(z0>-20)=NaN;
    vtop0(z0>-20)=NaN;
    utop(z0>-5)=NaN;
    vtop(z0>-5)=NaN;
    uu = interp2(LON0,LAT0,utop0,LON,LAT);
    vv = interp2(LON0,LAT0,vtop0,LON,LAT);
    uu2 = interp2(LON0,LAT0,utop,LON2,LAT2);
    vv2 = interp2(LON0,LAT0,vtop,LON2,LAT2);
    uu(LON > -124.25 & LAT > 47) = NaN;
    vv(LON > -124.25 & LAT > 47) = NaN;

	% For particular data sources
    if iopt<=3
		% and plot them
		subplot(131)
		set(gca,'Position',[-0.035 0.08 0.4 0.87])
		%    pcolor(LON0,LAT0,s0);
		%    shading interp
		gr_plot2(gr2.hgrid,sal,[slo shi]);
		caxis([slo shi]);
		colorbar('north')
		axis(aa); set(gca,'dataaspectratio',dar);
		hold on;
		% add lat/lon grid lines
		for ii=1:10
			plot([-126.00+0.25*ii -126.00+0.25*ii],[40 50],'w--','LineWidth',0.1);hold on; 
		end
		for ii=1:15
			plot([-130 -120],[45.0+0.25*ii 45.0+0.25*ii],'w--','LineWidth',0.1);hold on; 
		end
		% add velocity vectors
		ufact = .1;
		quiver(LON,LAT,ufact*uu*dar(1),ufact*vv,0,'k');
		% add a velocity scale
		uscale = 0.5;
		hold on; quiver(-123.75,45.65,ufact*uscale*dar(1),ufact*0,0,'k');
		hold on; quiver(-123.75,45.65,ufact*0*dar(1),ufact*uscale,0,'k');
		text(-123.9,45.6,[num2str(uscale),' m s^{-1}'],'color','k');
		% add salinity contours
		contour(LON0,LAT0,s0,sclevs,'-k');
		% add bathymetry contours
		[cc,hh] = contour(LON0,LAT0,z0,[-50 -50],'-w');
		for ii = 1:length(hh); set(hh(ii),'linewidth',2); end;
		[cc,hh] = contour(LON0,LAT0,z0,clevs,'-w');
		for ii = 1:length(hh); set(hh(ii),'linewidth',2); end;
		[cc,hh] = contour(LON0,LAT0,z0,[-1000 -1000],'-w');
		for ii = 1:length(hh); set(hh(ii),'linewidth',2); end;
		% add notes and labels
		title('(a) Surface Salinity (contours 30:0.5:34)');
		xlabel('Longitude (deg)');
		ylabel('Latitude (deg)');
		% add section line
	   	plot([slonw(1) slone(1)],[slat(1) slat(1)],'-k','linewidth',2)
	    plot([slonw(2) slone(2)],[slat(2) slat(2)],'-k','linewidth',2)
	    plot([slonw(3) slone(3)],[slat(3) slat(3)],'-k','linewidth',2)
	 	% add section text
	    text(-123.87,slat(nn),sect)
		old on;plot(x1,y1,'k-','LineWidth',2);hold on;grid on;
		lot([bb(1) bb(2)],[bb(3) bb(3)],'b-','LineWidth',2);hold on;
		lot([bb(1) bb(2)],[bb(4) bb(4)],'b-','LineWidth',2);hold on;
		lot([bb(1) bb(1)],[bb(3) bb(4)],'b-','LineWidth',2);hold on;
		lot([bb(2) bb(2)],[bb(3) bb(4)],'b-','LineWidth',2);hold on;
		ext(-125.48,45.1,'Isobaths = -1000  -200  -50 (m)','FontSize',14,'FontWeight','b','Color','k');hold on;
		ext(-125.48,45.2,ctime,'FontSize',14,'FontWeight','b','Color','k');hold on;
		ext(-125.48,45.3,['Forecast: ' modnames{imod}],'FontSize',14,'FontWeight','b','Color','k');hold on;
	    plot(x2,y2,'k+','LineWidth',3,'MarkerSize',14);hold on;
	    plot(xx2(nnn),yy2(nnn),'kO','LineWidth',3,'MarkerSize',14);hold on;
	    plot(x3,y3,'k-','LineWidth',3);hold on;
	    plot(x4,y4,'k-','LineWidth',3);hold on;

		subplot(132)
		set(gca,'Position',[0.3 0.08 0.25 0.87])
		%    pcolor(LON0,LAT0,s0);
		%    shading interp
		gr_plot2(gr2.hgrid,sal,[slo shi]);
		caxis([slo shi]);
		colorbar('north')
		axis(bb); set(gca,'dataaspectratio',dar);
		hold on;
		% add lat/lon grid lines
		for ii=1:5
			plot([-124.5+0.1*ii -124.5+0.1*ii],[40 50],'w--','LineWidth',0.1);hold on; 
		end
		for ii=1:7
			plot([-130 -120],[45.8+0.1*ii 45.8+0.1*ii],'w--','LineWidth',0.1);hold on; 
		end
		% add velocity vectors
		ufact = .1;
		quiver(LON2,LAT2,ufact*uu2*dar(1),ufact*vv2,0,'k');
		% add a velocity scale
		uscale = 0.5;
		hold on; quiver(-124.45,46.05,ufact*uscale*dar(1),ufact*0,0,'k');
		hold on; quiver(-124.45,46.05,ufact*0*dar(1),ufact*uscale,0,'k');
		text(-124.48,46.00,[num2str(uscale),' m s^{-1}'],'color','k');
		% add salinity contours
		contour(LON0,LAT0,s0,sclevs,'-k');
		% add bathymetry contours
		[cc,hh] = contour(LON0,LAT0,z0,[-50 -50],'-w');
		for ii = 1:length(hh); set(hh(ii),'linewidth',2); end;
		[cc,hh] = contour(LON0,LAT0,z0,clevs,'-w');
		for ii = 1:length(hh); set(hh(ii),'linewidth',2); end;
		[cc,hh] = contour(LON0,LAT0,z0,[-100 -100],'-w');
		for ii = 1:length(hh); set(hh(ii),'linewidth',2); end;
		% add notes and labels
		title('(b) Surface Salinity (contours 30:0.5:34)');
		xlabel('Longitude (deg)');
		ylabel('Latitude (deg)');
		% add section line
		%    plot([slonw(1) slone(1)],[slat(1) slat(1)],'-k','linewidth',2)
		%    plot([slonw(2) slone(2)],[slat(2) slat(2)],'-k','linewidth',2)
		%    plot([slonw(3) slone(3)],[slat(3) slat(3)],'-k','linewidth',2)
		% and section text
		%    text(-123.87,slat(nn),sect)
		hold on;plot(x1,y1,'k-','LineWidth',2);hold on;grid on;
		text(-124.48,45.83,'Isobaths = -200 -100  -50 (m)','FontSize',14,'FontWeight','b','Color','k');hold on;
		text(-124.48,45.86,ctime,'FontSize',14,'FontWeight','b','Color','k');hold on;
		text(-124.48,45.89,['Forecast: ' modnames{imod}],'FontSize',14,'FontWeight','b','Color','k');hold on;
		plot(x2,y2,'k+','LineWidth',3,'MarkerSize',14);hold on;
		%    plot(xx2(nnn),yy2(nnn),'kO','LineWidth',3,'MarkerSize',14);hold on;
		%    plot(x3,y3,'k-','LineWidth',3);hold on;
		%    plot(x4,y4,'k-','LineWidth',3);hold on;
	end
    
	if 0<1
    if iopt<=3
     subplot(233)
     set(gca,'Position',[0.6 0.5 0.37 0.52])
     title('(c) Surface Salinity in the Columbia River');
    else
     subplot(211)
     title('(a) Surface Salinity in the Columbia River');
    end
%    pcolor(LON0,LAT0,s0);
%    shading interp
    gr_plot2(gr2.hgrid,sal,[0 34]);
    caxis([0 34]);
    colorbar('north')
    axis(ee); set(gca,'dataaspectratio',dar);
    hold on;
    % add velocity vectors
%    ufact = .1;
%    quiver(LON2,LAT2,ufact*uu2*dar(1),ufact*vv2,0,'k');
    % add a velocity scale
%    uscale = 0.5;
%    hold on; quiver(-124.45,46.05,ufact*uscale*dar(1),ufact*0,0,'k');
%    hold on; quiver(-124.45,46.05,ufact*0*dar(1),ufact*uscale,0,'k');
%    text(-124.48,46.00,[num2str(uscale),' m s^{-1}'],'color','k');
    % add salinity contours
%    contour(LON0,LAT0,s0,sclevs2,'-k');
    % add bathymetry contours
    [cc,hh] = contour(LON0,LAT0,z0,[-10 -10],'-w');
    for ii = 1:length(hh); set(hh(ii),'linewidth',1.0); end;
%    [cc,hh] = contour(LON0,LAT0,z0,[-20 -20],'-w');
%    for ii = 1:length(hh); set(hh(ii),'linewidth',0.5); end;
    % add notes and labels
    xlabel('Longitude (deg)');
    ylabel('Latitude (deg)');
    hold on;plot(x1,y1,'k-','LineWidth',2);hold on;grid on;
    text(-123.75,46.15,['Forecast: ' modnames{imod}],'FontSize',14,'FontWeight','b','Color','k');hold on;
    text(-123.75,46.13,ctime,'FontSize',14,'FontWeight','b','Color','k');hold on;
    text(-123.75,46.11,'Isobaths = -10 (m)','FontSize',14,'FontWeight','b','Color','k');hold on;
%    plot(x2,y2,'k+','LineWidth',3,'MarkerSize',14);hold on;
    if iopt<=3
     subplot(236)
     set(gca,'Position',[0.6 0.01 0.37 0.52])
     title('(d) Bottom Salinity in the Columbia River');
    else
     subplot(212)
     title('(b) Bottom Salinity in the Columbia River');
    end
%    pcolor(LON0,LAT0,s0);
%    shading interp 
    gr_plot2(gr2.hgrid,salb,[0 34]);
    caxis([0 34]);
    colorbar('north')
    axis(ee); set(gca,'dataaspectratio',dar);
    hold on;
    % add velocity vectors
%    ufact = .1;
%    quiver(LON2,LAT2,ufact*uu2*dar(1),ufact*vv2,0,'k');
    % add a velocity scale
%    uscale = 0.5;
%    hold on; quiver(-124.45,46.05,ufact*uscale*dar(1),ufact*0,0,'k');
%    hold on; quiver(-124.45,46.05,ufact*0*dar(1),ufact*uscale,0,'k');
%    text(-124.48,46.00,[num2str(uscale),' m s^{-1}'],'color','k');
    % add salinity contours
%    contour(LON0,LAT0,s0,sclevs2,'-k');
    % add bathymetry contours
    [cc,hh] = contour(LON0,LAT0,z0,[-10 -10],'-w');
    for ii = 1:length(hh); set(hh(ii),'linewidth',1.0); end;
%    [cc,hh] = contour(LON0,LAT0,z0,[-20 -20],'-w');
%    for ii = 1:length(hh); set(hh(ii),'linewidth',0.5); end;
    % add notes and labels
    xlabel('Longitude (deg)');
    ylabel('Latitude (deg)');
    hold on;plot(x1,y1,'k-','LineWidth',2);hold on;grid on;
    text(-123.75,46.15,['Forecast: ' modnames{imod}],'FontSize',14,'FontWeight','b','Color','k');hold on;
    text(-123.75,46.13,ctime,'FontSize',14,'FontWeight','b','Color','k');hold on;
    text(-123.75,46.11,'Isobaths = -10 (m)','FontSize',14,'FontWeight','b','Color','k');hold on;
%    plot(x2,y2,'k+','LineWidth',3,'MarkerSize',14);hold on;
end

    % SELFE Salinity Section
    % set limits for the section plotting

  if 0>1
   for nn = 1:length(t_list);
    sect = sect_list{nn};
    switch sect
        case 'LP'
            ifg=333;fa='(c)';
        case 'GH'
            ifg=333;fa='(c)';
        case 'CR'
            ifg=336;fa='(d)';
        case 'CM'
            ifg=339;fa='(e)';
    end
    ob=ob_ini_fromTrasect(gr,t_list{nn});
    [lon,lat] = convm2ll(ob.xy.x,ob.xy.y);
    trLen = (slonw:(slone-slonw)/(length(ob.xy.x)-1):slone);
    ds=sz_readTimeStep(h,64);
    dsd=map_sz2hts(h, ds);
    sal_full=ob.xy.H*double(dsd);
    he=sz_readHeader([s_list{nn} 'elev.61']);
    de=sz_readTimeStep(he,it);
    e=ob.xy.H*double(de);
    slon = zeros(size(sal_full));
    for i = 1:size(sal_full,2);
        slon(:,i) = lon; %trLen';
    end
    dp=ob.xy.H*gr.hgrid.depth;
    sz_full=sz_computeZlevels(dp,e,gr.vgrid);
    z = -sum(gr.hgrid.depth(ob.xy.pidx)'.*ob.xy.w');
    % and plot it
    subplot(ifg);
    set(gca,'Position',scalea{nn});
    set(gca,'dataaspectratio',dar);
    pcolor(slon,sz_full,sal_full);
    caxis([slo shi]);
    axis([slonw(nn) slone(nn) -50 5]);
    shading interp;
    hold on
    % add salinity contours
    contour(slon,sz_full,sal_full,sclevs,'-k');
    % add bathymetry line
    plot(slon(:,1),z,'-k','linewidth',2)
    % add labels
    ylabel('Z (m)');
    title([fa,' ',sect,' Line Salinity']);
  end % loop: nn
 end
    frame=getframe(gcf);
    aviobj=addframe(aviobj,frame);
	%    print('-dpng', 'filetest.png');
 end % loop: it
end % loop: imod
aviobj=close(aviobj);
close all;
%!rm hgrid.gr3
%!rm hgrid.ll
%!rm ?_????.6?

