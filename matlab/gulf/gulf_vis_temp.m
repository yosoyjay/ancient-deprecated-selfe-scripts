close all
path(path,'./m-elio');
%--Coastline Information
load 'c17843.dat';
c1=c17843;
x1=c1(:,1);y1=c1(:,2);


start_day = datenum(2010, 1, 1, 0, 0, 0);
curr_day=datenum(yr, imn, idy, 0, 0, 0);
doy=curr_day-start_day;

% prepare fields - salt & veolocity
gr.hgrid=gr_readHGrid('hgrid.gr3');
gr2.hgrid=gr_readHGrid('hgrid.ll');
fg = hgrid2fg('hgrid.ll');

h=sz_readHeader('3_temp.63');
gr.vgrid=h.vgrid;
hv=sz_readHeader('3_hvel.64');

% set some default values
aa = [47.472609 63.592969 19.899843 30.46383]; % map region

r1pos = [46.659751359766 64.701418640234 19.899843 30.46383];
r2pos = [47.722045426667 55.634429426667 22.577343433332 30.489727433332];
r3pos = [49.832014493333 57.744398493333 21.601482739998 29.513866739998];
r4pos = [55 58 25 27.5];

%fn = sprintf('gulf-sal-%d-%02d-%02d.avi', 2010, 12, 20);
%aviobj=avifile(fn,'fps',0.5);

ist=1;intv=1;ied=96;
vbot = 18;

for it=ist:intv:ied 
ctime = [datestr(curr_day + (it * 900.0 / 86400.0), 'yyyy-mm-dd HH:MM:SS') ' UTC'];
%it = 1
[d ts]=sz_readTimeStep(h,it); %step 64 is 1600 hrs PST (= 0000 GMT)
% prepare fields - velocity
[dh tsh]=sz_readTimeStep(hv,it); %step 64 is 1600 hrs PST (= 0000 GMT)
u = map_sz2hts(hv,dh(:,1),1);
v = map_sz2hts(hv,dh(:,2),1);
u = u(:,end);
v = v(:,end);

[s] = map_sz2hts(h,d,1);
[s0] = s(:,end); % End is the surface layer
[s1] = s(:,vbot); % vbot is the bottom layer
s0(s0<0) = NaN;
s1(s1<0) = NaN;
sal=s0;
salb=s1;

slo = 10; shi = 30;

myfig = figure;
figure(myfig)
set(gcf,'position',[0 0 1280 720],'Color','w');
fs1 = 14; % fontsize
set(gcf,'defaultaxesfontsize',fs1);
set(gcf,'defaulttextfontsize',fs1);
set(gcf,'defaulttextfontweight','b');
set(gcf,'defaultaxesfontweight','b');
set(gcf,'defaultaxeslinewidth',2);
set(gcf,'defaultaxesfontname','Times');
set(gcf,'defaulttextfontname','Times');

ufact = 1.0;
uscale = 1.0;

% region 1
ax1 = axes;
set(ax1,'dataaspectratio',[1 1 1]);
set(ax1,'Position',[0.1 0.1 0.85 0.9])
axis(r1pos);
% 46.659751359766+19.899843+64.701418640234+30.46383
caxis([slo shi]);
cb = colorbar('EastOutside')
grid(ax1, 'on');
box(ax1, 'on');
title(['Surface Temperature at ' ctime]);
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
gr_plot2(gr2.hgrid,sal,[slo shi]);
hold on;
plot(x1,y1,'k-','LineWidth',2);hold on;
%quiver(hv.hgrid.x(1:10:end),hv.hgrid.y(1:10:end),u(1:10:end),v(1:10:end),0,'k');
hold on;
%lon=54.79964, lat=22.48802 ,
%quiver(48,22,ufact*uscale,ufact*0,0,'k');hold on;
%quiver(54.79964, 22.48802,ufact*0,ufact*uscale,0,'k');hold on;
%text(48, 21.5,[num2str(uscale),' m s^{-1}'],'color','k');hold on;

%ax2 = axes;
%grid(ax2, 'on');
%box(ax2, 'on');
%set(ax2,'dataaspectratio',[1 1 1]);
%set(ax2,'Position',[0.75 0.7 0.3 0.3])
%axis([47.722045426667 55.634429426667 22.577343433332 30.489727433332]);
%axis(r2pos);
%%47.722045426667+22.577343433332+55.634429426667+30.489727433332
%caxis([slo shi]);
%gr_plot2(gr2.hgrid,sal,[slo shi]);
%hold on;
%plot(x1,y1,'k-','LineWidth',2);hold on;
%%quiver(hv.hgrid.x(1:5:end),hv.hgrid.y(1:5:end),u(1:5:end),v(1:5:end),0,'k');
%hold on;

%ax3 = axes;
%grid(ax3, 'on');
%box(ax3, 'on');
%set(ax3,'dataaspectratio',[1 1 1]);
%set(ax3,'Position',[0.75 0.4 0.3 0.3])
%axis(r3pos);
%%47.722045426667+22.577343433332+55.634429426667+30.489727433332
%caxis([slo shi]);
%gr_plot2(gr2.hgrid,sal,[slo shi]);
%hold on;
%plot(x1,y1,'k-','LineWidth',2);hold on;
%%quiver(hv.hgrid.x(1:5:end),hv.hgrid.y(1:5:end),u(1:5:end),v(1:5:end),0,'k');
%hold on;

ax4 = axes;
grid(ax4, 'on');
box(ax4, 'on');
set(ax4,'dataaspectratio',[1 1 1]);
set(ax4,'Position',[0.6 0.6 0.25 0.25])
set(ax4,'FontSize',10);
axis(r4pos);
%47.722045426667+22.577343433332+55.634429426667+30.489727433332
caxis([slo shi]);
gr_plot2(gr2.hgrid,sal,[slo shi]);
hold on;
plot(x1,y1,'k-','LineWidth',2);hold on;
quiver(hv.hgrid.x,hv.hgrid.y,u,v,0,'k');
%quiver(hv.hgrid.x(1:5:end),hv.hgrid.y(1:5:end),u(1:5:end),v(1:5:end),0,'k');
hold on;
fn = sprintf('anim/gulf-tem-%d-%02d-%02d-%03d.png', 2010, 12, 20, it);
print('-dpng', fn);
%frame=getframe(gcf);
%aviobj=addframe(aviobj,frame);
close(myfig);

end % it
