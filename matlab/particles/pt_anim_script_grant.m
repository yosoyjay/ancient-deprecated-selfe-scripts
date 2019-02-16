pathset = '_swim';
rec_vid = 'y';% Yes ('y'), or no ('n')
xlim = [-124.0905 -123.6661];
ylim = [46.1317   46.3186];

fg = hgrid2fg('hgrid.gr3');
fgll = fg;
[fgll.x,fgll.y] = convm2ll(fgll.x,fgll.y);

%path = partpath(['particle' pathset '.pth']);
path=partpath('particle.pth');
[path.lon,path.lat] = convm2ll(path.x,path.y);

[parts,steps] = size(path.lon);
%steps = steps-1;

fid = fopen('dry_elements.dat');
fscanf(fid,'%s',7);
for i = 1:steps
    count = fscanf(fid,'%i',1);
    eval(['drydata.step' num2str(i) ' = fscanf(fid,''%i\n'',count);'])
end

figure
set(gcf,'Position',[55 158 1159 707])
set(gca,'XLim',xlim)
set(gca,'YLim',ylim)

hold on
elems = drawelems(fgll);
bnd = plotbnd(fgll);
dryspots = patch(fgll.x(fgll.e(1,:)),fgll.y(fgll.e(1,:)),[-.1 -.1 -.1],[.5 .5 .5]);
set(dryspots,'EdgeColor','none')
tracers = 1;
tracer.x = ones(parts,steps,tracers);
tracer.y = ones(parts,steps,tracers);
for i = 1:tracers
    tracer.x(:,:,i) = [ones(parts,i).*nan path.lon(:,1:steps-i)];
    tracer.y(:,:,i) = [ones(parts,i).*nan path.lat(:,1:steps-i)];
end

for i = tracers:-1:1
    cdata = [i/(2*tracers) i/(2*tracers) i/(2*tracers)];
    %sdata = 6 * (1-i/(tracers+1));
    sdata = 6;
    eval(['part' num2str(i) ' = plot(squeeze(tracer.x(1,1,' num2str(i) ')),squeeze(tracer.y(1,1,' num2str(i) ')),''LineStyle'',''.'',''Marker'',''.'',''Color'',cdata,''MarkerSize'',sdata);'])
end
particles = scatter(path.lon(:,1),path.lat(:,1),12,[zeros(parts,1) zeros(parts,1) zeros(parts,1)],'filled');

if strcmp(rec_vid,'y')
    mov = avifile('part_path.avi');
end
flag = 1;
while flag==1
   for i = 1:steps
        eval(['dry = drydata.step' num2str(i) ';']);
        temp1 = fgll.e(dry,:);
        temp2 = fgll.x(temp1);
        x_dat = fgll.x(fgll.e(dry,:));
        y_dat = fgll.y(fgll.e(dry,:));
        z_dat = ones(3,length(dry));
        set(dryspots,'XData',fgll.x(fgll.e(dry,:))','YData',fgll.y(fgll.e(dry,:))','ZData',ones(3,length(dry)));
        zcol = (1-exp(path.z(:,i))).*0.8;
        set(particles,'XData',path.lon(:,i),'YData',path.lat(:,i),'CData',[zcol zcol zcol],'SizeData',(1-zcol).*12)
        for j = 1:tracers
            eval(['set(part' num2str(j) ',''XData'',squeeze(tracer.x(:,i,j)),''YData'',squeeze(tracer.y(:,i,j)));'])
        end
        drawnow
        if strcmp(rec_vid,'y')
            F = getframe(gca);
            mov = addframe(mov,F);
            if i==steps
                flag = 0;
            end
        end
    end
end
if strcmp(rec_vid,'y')
    mov = close(mov);
end
