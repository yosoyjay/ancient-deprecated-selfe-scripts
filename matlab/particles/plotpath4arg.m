function h = plotpath(format,aux,titlestr,clustsz)

pathdata = partpath('particle.pth');
fg = hgrid2fg('hgrid.gr3');
fgll = hgrid2fg('hgrid.ll');
fgll.z = fg.z;
fgll.bnd = fg.bnd;
fgll.bndpth = fg.bndpth;

particles = size(pathdata.x,1);
% Convert particle positions to lat/lon
pathdata.lon = zeros(size(pathdata.x));
pathdata.lat = pathdata.lon;
for i = 1:size(pathdata.x,1)
    for j = 1:size(pathdata.x,2)
        dist = sqrt( (fg.x-pathdata.x(i,j)).^2 + (fg.y-pathdata.y(i,j)).^2 );
        [d,ind] = sort(dist);
        d = d(1:4);
        ind = ind(1:4);
        xl = min(fg.x(ind));
        xli = find(fg.x==xl);xli=xli(1);
        xh = max(fg.x(ind));
        xhi = find(fg.x==xh);xhi=xhi(1);
        xrat = (pathdata.x(i,j)-fg.x(xli)) / (fg.x(xhi)-fg.x(xli));
        yl = min(fg.y(ind));
        yli = find(fg.y==yl);yli = yli(1);
        yh = max(fg.y(ind));
        yhi = find(fg.y==yh);yhi=yhi(1);
        yrat = (pathdata.y(i,j)-fg.y(yli)) / (fg.y(yhi)-fg.y(yli));
        pathdata.lon(i,j) = fgll.x(xli) + xrat*(fgll.x(xhi)-fgll.x(xli));
        pathdata.lat(i,j) = fgll.y(yli) + yrat*(fgll.y(yhi)-fgll.y(yli));
    end
end

xmin = min(min(pathdata.lon))-.01;
xmax = max(max(pathdata.lon))+.01;
ymin = min(min(pathdata.lat))-.01;
ymax = max(max(pathdata.lat))+.01;
xrange = xmax-xmin;
yrange = ymax-ymin;
if xrange>yrange
    rat = yrange/xrange;
    xpix = 1100;
    ypix = max([1100*rat 800]);
    rngadd = (max([xrange*8/11 yrange])-yrange)/2;
    ymin = ymin - rngadd;
    ymax = ymax + rngadd;
else
    rat = xrange/yrange;
    xpix = max([1100*rat 800]);
    ypix = 1100;
    rngadd = (max([yrange*8/11 xrange])-xrange)/2;
    xmin = xmin - rngadd;
    xmax = xmax + rngadd;
end
moving = find(pathdata.x(1,:)~=pathdata.x(1,1));
reltime = pathdata.time(moving(1)-1);
startstr = [num2str(floor(reltime/3600),'%2.2i') ':' num2str(rem(reltime/3600,1)*60,'%2.2i')];

switch true
    case strcmp(format(1:4),'simu')|strcmp(format(1:4),'Simu')
        h.fig = figure('Position',[20 20 xpix ypix]);
        set(gca,'Color',[.2 .2 .2])
        hold on
        
        seafloor = patch('Vertices',[fgll.x fgll.y fgll.z],'Faces',fgll.e,'FaceVertexCData',...
            ones(length(fgll.x),3).*1,'FaceColor','interp','EdgeColor','none',...
            'SpecularStrength',0,'DiffuseStrength',0.7);
        L = light;
        lighting phong
        lightangle(L,180,60)
        material dull
        h.bnd = plotbnd(fgll);
        h.path = plot(pathdata.lon',pathdata.lat','r');
        h.min15 = plot(pathdata.lon,pathdata.lat,'.r');
        h.stpnt = plot(pathdata.lon(:,1),pathdata.lat(:,1),'ok');
        h.endpnt = plot(pathdata.lon(:,96),pathdata.lat(:,96),'.k');
        h.hourly = plot(pathdata.lon(:,4:4:96),pathdata.lat(:,4:4:96),'ok');
        east = find(pathdata.lon(:,1)==max(pathdata.lon(:,1)));
        h.startstr = text(pathdata.lon(east)+xrange/60,pathdata.lat(east),num2str(startstr));
        set(gca,'XLim',[xmin xmax],'YLim',[ymin ymax])
        set(gca,'ZLim',[-100 100]);
        set(gca,'FontSize',12)
        xlabel('Longitude')
        ylabel('Latitude')
        title(titlestr)
        box on
        % Draw gridlines
        xg = [get(gca,'XTick');get(gca,'XTick')];
        yg = [get(gca,'YTick');get(gca,'YTick')];
        xg2 = [xg [ones(1,length(yg)).*-130;ones(1,length(yg)).*-100]];
        yg2 = [[ones(1,length(xg)).*30;ones(1,length(xg)).*60] yg];
        plot3(xg2,yg2,ones(size(xg2)).*20,':k');
        
        % Write trajectories to a text file
        hour = floor(pathdata.time(moving(1)-1:end)./3600);
        minute = (pathdata.time(moving(1)-1:end) - hour.*3600)./60;
        latdeg = floor(pathdata.lat(moving(1)-1:end));
        latdec = mod(pathdata.lat(moving(1)-1:end),1).*60;
        londeg = floor(pathdata.lon(moving(1)-1:end));
        londec = mod(pathdata.lon(moving(1)-1:end),1).*60;
        textdata = [hour' minute' latdeg' latdec' londeg' londec'];
        fid = fopen('drifterpath.dat','w');
        fprintf(fid,'%2.2i%2.2i   %2.2i %6.4f %2.2i %6.4f\n',textdata');
        fclose(fid);
    case strcmp(format(1:4),'tran')|strcmp(format(1:4),'Tran')
        h.fig = figure('Position',[20 20 xpix ypix]);
        set(gca,'Color',[.2 .2 .2])
        hold on
        
        seafloor = patch('Vertices',[fgll.x fgll.y fgll.z],'Faces',fgll.e,'FaceVertexCData',...
            ones(length(fgll.x),3).*1,'FaceColor','interp','EdgeColor','none',...
            'SpecularStrength',0,'DiffuseStrength',0.7,'FaceAlpha',0.5);
        L = light;
        lighting phong
        lightangle(L,180,60)
        material dull
        plotbnd(fgll);
        plot(pathdata.lon',pathdata.lat','r');
        plot(pathdata.lon,pathdata.lat,'.r');
        plot(pathdata.lon(:,1),pathdata.lat(:,1),'ok');
        plot(pathdata.lon(:,96),pathdata.lat(:,96),'.k');
        plot(pathdata.lon(:,96),pathdata.lat(:,96),'.k')
        plot(pathdata.lon(:,4:4:96),pathdata.lat(:,4:4:96),'ok')
        east = find(pathdata.lon(:,1)==max(pathdata.lon(:,1)));
        text(pathdata.lon(east)+xrange/60,pathdata.lat(east),num2str(startstr));
        grid
        grid minor
        set(gca,'XLim',[xmin xmax],'YLim',[ymin ymax])
        set(gca,'ZLim',[-100 100]);
        set(gca,'FontSize',12)
        xlabel('Longitude')
        ylabel('Latitude')
        title(titlestr)
        box on
        % Draw gridlines
        xg = [get(gca,'XTick');get(gca,'XTick')];
        yg = [get(gca,'YTick');get(gca,'YTick')];
        xg2 = [xg [ones(1,length(yg)).*-130;ones(1,length(yg)).*-100]];
        yg2 = [[ones(1,length(xg)).*30;ones(1,length(xg)).*60] yg];
        plot3(xg2,yg2,ones(size(xg2)).*20,':k');
        
    case strcmp(format(1:4),'comp')|strcmp(format(1:4),'Comp')
        % Figure out how many releases are present
        reps = sum(pathdata.x(:,1)==pathdata.x(1,1)&pathdata.y(:,1)==pathdata.y(1,1));
        releases = round(size(pathdata.x,1)./reps);
        relpos.x = aux(:,1); relpos.y = aux(:,2);
        relpos.lon = zeros(size(relpos.x));
        relpos.lat = relpos.lon;
        for i = 1:size(relpos.x,1)
            for j = 1:size(relpos.x,2)
                dist = sqrt( (fg.x-relpos.x(i,j)).^2 + (fg.y-relpos.y(i,j)).^2 );
                [d,ind] = sort(dist);
                d = d(1:4);
                ind = ind(1:4);
                xl = min(fg.x(ind));
                xli = find(fg.x==xl);xli=xli(1);
                xh = max(fg.x(ind));
                xhi = find(fg.x==xh);xhi=xhi(1);
                xrat = (pathdata.x(i,j)-fg.x(xli)) / (fg.x(xhi)-fg.x(xli));
                yl = min(fg.y(ind));
                yli = find(fg.y==yl);yli = yli(1);
                yh = max(fg.y(ind));
                yhi = find(fg.y==yh);yhi=yhi(1);
                yrat = (relpos.y(i,j)-fg.y(yli)) / (fg.y(yhi)-fg.y(yli));
                relpos.lon(i,j) = fgll.x(xli) + xrat*(fgll.x(xhi)-fgll.x(xli));
                relpos.lat(i,j) = fgll.y(yli) + yrat*(fgll.y(yhi)-fgll.y(yli));
            end
        end
        drops = size(relpos.x,1);
        brad = clustsz.big;% Includes all particles released within clustsz.big meters
        srad = clustsz.small;% Includes all particles released within clustsz.small meters
        bpnts = ceil(brad/10);
        spnts = ceil(srad/10);
        brho=ones(1,bpnts)*brad;
        srho=ones(1,spnts)*srad;
        btheta=linspace(0,2*pi,bpnts);
        stheta=linspace(0,2*pi,spnts);
        [bcirx,bciry] = pol2cart(btheta,brho);
        [scirx,sciry] = pol2cart(stheta,srho);
        
        % Determine the times of the releases
        for i = 1:reps
            moving = find(pathdata.x((i-1)*releases+1,:)~=pathdata.x((i-1)*releases+1,1));
            reltime(i) = pathdata.time(moving(1)-1);
        end

        h = figure('Position',[20 20 xpix ypix]);
        xmin = min(min(pathdata.lon))-.1;
        xmax = max(max(pathdata.lon))+.1;
        ymin = min(min(pathdata.lat))-.1;
        ymax = max(max(pathdata.lat))+.1;
        cnt = 0;
        for i = 1:drops
            bpolx=bcirx+relpos.x(i);
            bpoly=bciry+relpos.y(i);
            spolx=scirx+relpos.x(i);
            spoly=sciry+relpos.y(i);
            inbig = inpolygon(pathdata.x(1:releases,1),pathdata.y(1:releases,2),bpolx,bpoly);
            insml = inpolygon(pathdata.x(1:releases,1),pathdata.y(1:releases,2),spolx,spoly);
            for j = 1:reps
                cnt = cnt+1;
                tmpx = pathdata.lon((j-1)*releases+1:j*releases,:);
                tmpy = pathdata.lat((j-1)*releases+1:j*releases,:);
                subplot(drops,reps,cnt)
                set(gca,'Color',[.2 .2 .2])
                set(gca,'FontSize',12)
                hold on
                set(gca,'XTick',[],'YTick',[])
                
                seafloor = patch('Vertices',[fgll.x fgll.y fgll.z],'Faces',fgll.e,'FaceVertexCData',...
                    ones(length(fgll.x),3).*1,'FaceColor','interp','EdgeColor','none',...
                    'SpecularStrength',0,'DiffuseStrength',0.7,'FaceAlpha',.5);
             	L = light;
                lighting phong
                lightangle(L,180,60)
             	material dull
                plotbnd(fgll);
                plot(tmpx(inbig,:)',tmpy(inbig,:)','r')
                plot(tmpx(insml,:)',tmpy(insml,:)','b')
                set(gca,'XLim',[xmin xmax]);
                set(gca,'YLim',[ymin ymax]);
                set(gca,'ZLim',[-100 100]);
                box on
                if j==1
                    labstr = {[num2str(floor(relpos.lat(i))) ' ' num2str(rem(relpos.lat(i),1)*60,'%6.4f')];[num2str(floor(relpos.lon(i))) ' ' num2str(-rem(relpos.lon(i),1)*60,'%6.4f')]};
                    ylabel(labstr)
                end
                if i==drops
                    xlabel([num2str(floor(reltime(j)/3600),'%2.2i') ':' num2str(rem(reltime(j)/3600,1)*60,'%2.2i')]);
                end
            end
        end
        subplot(drops,reps,ceil(reps/2))
        title(titlestr);
end
