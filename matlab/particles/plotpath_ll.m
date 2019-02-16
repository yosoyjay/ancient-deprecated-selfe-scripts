function h = plotpath(format,aux,reltime,titlestr,clustsz,ipt,obs,obs_times)

pathdata = partpath_ll('particle.ll');
fgll = hgrid2fg('hgrid.ll');
fg = hgrid2fg('hgrid.gr3');
bnd = input_complex_bnd('f22.bnd');

particles = size(pathdata.x,1);
% Convert particle positions to lat/lon
% [pathdata.x,pathdata.lat] = convm2ll(pathdata.x,pathdata.y);

if size(ipt, 1) == 0
    xmin = min(min(pathdata.x))-.01;
    xmax = max(max(pathdata.x))+.01;
    ymin = min(min(pathdata.y))-.01;
    ymax = max(max(pathdata.y))+.01;
else  % make dims more of a suggestion than a requirement
    xmin = min([min(min(pathdata.x)), ipt(1)])-0.01;
    xmax = max([max(max(pathdata.x)), ipt(1)])+0.01;
    ymin = min([min(min(pathdata.y)), ipt(2)])-0.01;
    ymax = max([max(max(pathdata.y)), ipt(2)])+0.01;
end
if (length(obs)>0)
    if (min(obs(:,1))<xmin)
        xmin = min(obs(:,1));
    end
    if (max(obs(:,1))>xmax)
        xmax = max(obs(:,1));
    end    
    if (min(obs(:,2))<ymin)
        ymin = min(obs(:,2));
    end
    if (max(obs(:,2))>ymax)
        ymax = max(obs(:,2));
    end       
end

xrange = xmax-xmin;
yrange = ymax-ymin;

long_side = 800;
short_side = 600;

if xrange>yrange
    rat = yrange/xrange;
    xpix = long_side;
    ypix = max([long_side*rat short_side]);
    rngadd = (max([xrange*short_side/long_side yrange])-yrange)/2;
    ymin = ymin - rngadd;
    ymax = ymax + rngadd;
else
    rat = xrange/yrange;
    xpix = max([long_side*rat short_side]);
    ypix = long_side;
    rngadd = (max([yrange*short_side/long_side xrange])-xrange)/2;
    xmin = xmin - rngadd;
    xmax = xmax + rngadd;
end

moving = find(pathdata.x(1,:)~=pathdata.x(1,1));
%reltime = pathdata.time(moving(1)-1);

switch true
    case strcmp(format(1:4),'simu')|strcmp(format(1:4),'Simu')
        startstr = [num2str(floor(reltime./3600),'%2.2d') ':' num2str(round(rem(reltime./3600,1)*60),'%2.2d')];
        h.fig = figure('Position',[20 20 xpix ypix]);
        set(gca,'Color',[1, 1, 1]);
        hold;
        plot(bnd(:,1), bnd(:,2), 'k');
%         colormap(gray);
%         set(gca,'CLim', [-200, 0]);
%         seafloor = patch('Vertices',[fgll.x fgll.y fgll.z],'Faces',fgll.e,'FaceVertexCData',...
%             fg.z,'FaceColor','interp','EdgeColor','none',...
%             'SpecularStrength',0,'DiffuseStrength',0.7,'FaceAlpha',0.5);

%         seafloor = patch('Vertices',[fgll.x fgll.y fgll.z],'Faces',fgll.e,'FaceVertexCData',...
%             ones(length(fgll.x),3).*1,'FaceColor','interp','EdgeColor','none',...
%             'SpecularStrength',0,'DiffuseStrength',0.7);
%         L = light;
%         lighting phong
%         lightangle(L,180,60)
%         material dull


%         h.bnd = plotbnd(fgll);
        h.path = plot(pathdata.x',pathdata.y','r');
        h.min15 = plot(pathdata.x,pathdata.y,'.r');
        h.stpnt = plot(pathdata.x(:,1),pathdata.y(:,1),'ok');
        h.endpnt = plot(pathdata.x(:,end),pathdata.y(:,end),'.k');
        firsthour = 4 - floor(mod(rem(reltime/3600,1)*4+1,4));
%         h.hourly = plot(pathdata.x(:,firsthour:4:end),pathdata.y(:,firsthour:4:end),'ok');
        h.hourly = scatter(pathdata.x(:,firsthour:4:end),pathdata.y(:,firsthour:4:end), 36,...
                           'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
                       
%         h.hourly = plot(pathdata.x(:,firsthour:4:end),pathdata.y(:,firsthour:4:end),'ok');
                       
        east = find(pathdata.x(:,1)==max(pathdata.x(:,1)));
        h.startstr = text(pathdata.x(east)+xrange/60,pathdata.y(east),num2str(startstr));
        set(gca,'XLim',[xmin xmax],'YLim',[ymin ymax])
%         set(gca,'ZLim',[-100 100]);
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
        % Check if there are any drifter data for this time period
        %['!wget -o t.log -O t.csv ''http://data.stccmop.org/ws/product/factory.py/getcsv?&endtime=' year '-' month '-' day 'T'];
        
        if length(obs)>0
            % fill in observations
%             otimes = (obs_times(1):15*60/86400:obs_times(end));
%             olon = interp(obs_times, obs(:,1), otimes);
%             olat = interp(obs_times, obs(:,2), otimes);
%             plot(olon, olat, 'g');
%             scatter(olon, olat, 'Marker', 'd', 'MarkerFaceColor', 'g');
%             firsthour = 
            plot(obs(:,1), obs(:,2), 'g');
%             scatter(obs(:,1), obs(:,2), 'Marker', 'd', 'MarkerFaceColor', 'g'); %, 'SizeData', 25);
            scatter(obs(:,1),obs(:,2), 36, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'b');
%             for tt = 1:2:length(obs_times)
%                 dvec = datevec(obs_times(tt));
%                 text(obs(tt,1), obs(tt,2), datestr(obs_times(tt), 'mm-dd-yy HH:MM'), 'Color', 'g');
%             end
        end
        
        %%%%%%%% Write trajectories to a text file
%         hour = floor(pathdata.time(moving(1)-1:end)./3600);
%         minute = (pathdata.time(moving(1)-1:end) - hour.*3600)./60;
%         latdeg = floor(pathdata.y(moving(1)-1:end));
%         latdec = rem(pathdata.y(moving(1)-1:end),1).*60;
%         londeg = ceil(pathdata.x(moving(1)-1:end));
%         londec = rem(-pathdata.x(moving(1)-1:end),1).*60;
%         textdata = [hour' minute' latdeg' latdec' londeg' londec'];
%         fid = fopen('drifterpath.dat','w');
%         fprintf(fid,'%2.2i%2.2i   %2.2i %6.3f %2.2i %6.3f\n',textdata');
%         fclose(fid);
    case strcmp(format(1:4),'tran')|strcmp(format(1:4),'Tran')
        h.fig = figure('Position',[20 20 xpix ypix]);
        set(gca,'Color',[1, 1, 1])
        hold;
        plot(bnd(:,1), bnd(:,2), 'k');
%         seafloor = patch('Vertices',[fgll.x fgll.y fgll.z],'Faces',fgll.e,'FaceVertexCData',...
%             ones(length(fgll.x),3).*1,'FaceColor','interp','EdgeColor','none',...
%             'SpecularStrength',0,'DiffuseStrength',0.7,'FaceAlpha',0.5);
%         L = light;
%         lighting phong
%         lightangle(L,180,60)
%         material dull
%         plotbnd(fgll);
        plot(pathdata.x',pathdata.y','r');
        plot(pathdata.x,pathdata.y,'.r');
        plot(pathdata.x(:,1),pathdata.y(:,1),'ok');
        plot(pathdata.x(:,end),pathdata.y(:,end),'.k');
        plot(pathdata.x(:,end),pathdata.y(:,end),'.k')
        plot(pathdata.x(:,4:4:end),pathdata.y(:,4:4:end),'ok')
        east = find(pathdata.x(:,1)==max(pathdata.x(:,1)));
        text(pathdata.x(east)+xrange/60,pathdata.y(east),num2str(startstr));
        grid
        grid minor
        set(gca,'XLim',[xmin xmax],'YLim',[ymin ymax])
        %set(gca,'ZLim',[-100 100]);
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
%         relpos.x = zeros(size(relpos.x));
%         relpos.y = relpos.x;
%         for i = 1:size(relpos.x,1)
%             for j = 1:size(relpos.x,2)
%                 dist = sqrt( (fg.x-relpos.x(i,j)).^2 + (fg.y-relpos.y(i,j)).^2 );
%                 [d,ind] = sort(dist);
%                 d = d(1:4);
%                 ind = ind(1:4);
%                 xl = min(fg.x(ind));
%                 xli = find(fg.x==xl);xli=xli(1);
%                 xh = max(fg.x(ind));
%                 xhi = find(fg.x==xh);xhi=xhi(1);
%                 xrat = (pathdata.x(i,j)-fg.x(xli)) / (fg.x(xhi)-fg.x(xli));
%                 yl = min(fg.y(ind));
%                 yli = find(fg.y==yl);yli = yli(1);
%                 yh = max(fg.y(ind));
%                 yhi = find(fg.y==yh);yhi=yhi(1);
%                 yrat = (relpos.y(i,j)-fg.y(yli)) / (fg.y(yhi)-fg.y(yli));
%                 relpos.x(i,j) = fgll.x(xli) + xrat*(fgll.x(xhi)-fgll.x(xli));
%                 relpos.y(i,j) = fgll.y(yli) + yrat*(fgll.y(yhi)-fgll.y(yli));
%             end
%         end
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
        xmin = min(min(pathdata.x))-.1;
        xmax = max(max(pathdata.x))+.1;
        ymin = min(min(pathdata.y))-.1;
        ymax = max(max(pathdata.y))+.1;
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
                tmpx = pathdata.x((j-1)*releases+1:j*releases,:);
                tmpy = pathdata.y((j-1)*releases+1:j*releases,:);
                subplot(drops,reps,cnt)
                set(gca,'Color',[.2 .2 .2])
                set(gca,'FontSize',12)
                hold on
                set(gca,'XTick',[],'YTick',[])

                plot(bnd(:,1), bnd(:,2), 'k');
%                 seafloor = patch('Vertices',[fgll.x fgll.y fgll.z],'Faces',fgll.e,'FaceVertexCData',...
%                     ones(length(fgll.x),3).*1,'FaceColor','interp','EdgeColor','none',...
%                     'SpecularStrength',0,'DiffuseStrength',0.7,'FaceAlpha',.5);
%              	L = light;
%                 lighting phong
%                 lightangle(L,180,60)
%              	material dull
%                 plotbnd(fgll);
                plot(tmpx(inbig,:)',tmpy(inbig,:)','r')
                plot(tmpx(insml,:)',tmpy(insml,:)','b')
                set(gca,'XLim',[xmin xmax]);
                set(gca,'YLim',[ymin ymax]);
                %set(gca,'ZLim',[-100 100]);
                box on
                if j==1
                    labstr = {[num2str(floor(relpos.y(i))) ' ' num2str(rem(relpos.y(i),1)*60,'%6.4f')];[num2str(ceil(relpos.x(i))) ' ' num2str(-rem(relpos.x(i),1)*60,'%6.4f')]};
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
