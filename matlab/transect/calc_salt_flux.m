function [var1outflow, var1inflow] = calc_salt_flux(runPath, runName)
% function [var1outflow, var1inflow] = calc_salt_flux(runPath, runName)
%  Calculates and returns the flux across a boundary.  Hacked up to use for salt flux across
% the mouth.
% runpath should be in a cell array.
% runname should be in a cell array too
 
path(path,'/usr/local/cmop/matlab/cmop/m-elio');
% Build points file
%basedir = '/home/workspace/ccalmr/hindcasts/2009-30-22/run/';
%basedir = '/home/workspace/local0/forecasts/f22/yesterday/run/';
%basedir_all = {'/home/workspace/project/lopezj/density/2011/baseline/run/',...
%               '/home/workspace/project/lopezj/density/2011/combo_b/run/'}

basedir_all = runPath;
runname_all = runName;

varnames = {'salt'};
var1name ='salt';
cax = [0 100];
ylims = [-20 3];
%runname_all = {'Baseline';'ComboB'};
endnum = [14; 21];
for ii = 1 %:length(basedir_all)
    first = 0;
    basedir = basedir_all{ii};
    runname = runname_all{ii};
    for dd = 5:8
        disp(['day',num2str(dd)]);
        h1=sz_readHeader([basedir,'/outputs/',num2str(dd),'_',varnames{1},'.63']);
%        h2=sz_readHeader([basedir,'outputs/',num2str(dd),'_',varnames{2},'.63']);
        hv=sz_readHeader([basedir,'/outputs/',num2str(dd),'_hvel.64']);
        he=sz_readHeader([basedir,'/outputs/',num2str(dd),'_elev.61']);
        if first == 0
            first = 1;
            gr.hgrid=gr_readHGrid([basedir,'/hgrid.gr3']);
            gr.vgrid=h1.vgrid;
            gr.hgrid.tri = gr.hgrid.elem(:,3:5);
            % compute transect
        %    [ob]= ob_ini_fromTrasect(gr, '/home/workspace/project/lopezj/density/2011/plots/db22_outer_mouth.bp');
            [ob]= ob_ini_fromTrasect(gr, '/home/workspace/project/lopezj/density/2011/plots/mouth_4.bp');
            badnodes = isnan(ob.xy.w(:,1));
            last = find(~isnan(ob.xy.w(:,1)), 1, 'last' );
            trDist = sqrt((ob.xy.x(2:end)-ob.xy.x(1:end-1)).^2+(ob.xy.y(2:end)-ob.xy.y(1:end-1)).^2);
            trLen = cumsum(trDist);
			if last > size(trLen,1)
				last = size(trLen,1);
			end
            
            x_comp = repmat( diff(ob.xy.y)./trDist , 1 , h1.vgrid.nLevels-1);
            y_comp = repmat( -diff(ob.xy.x)./trDist,  1 , h1.vgrid.nLevels-1);
            
            trDist = repmat(trDist, 1 , h1.vgrid.nLevels-1);
            trLen = repmat(trLen, 1 , h1.vgrid.nLevels-1);
            
            cross = [];
            along = [];
            
            inflow = [];
            outflow = [];
            
            var1inflow = [];
            var1outflow = [];
%            var2inflow = [];
%           var2outflow = [];
            
            t = [];
            i = 0;
        end;
        for step=1:90
%            if dd ==17 && step <8
%                continue
%            end
            if ~isempty(t)
                t(end+1) = t(end)+960/86400;
            else
                t = 960/86400;
            end
            i = i+1;
            % Read timestep of data
            %variable 1
            ds=sz_readTimeStep(h1,step);
            dsd=map_sz2hts(h1, ds);
            dry = dsd==-99;
            dsd(dry) = nan;
            
            var1=ob.xy.H*double(dsd);
            var1(badnodes,:) = nan;
%            % variable 2
            
            dv=sz_readTimeStep(hv,step);
            dud=map_sz2hts(h1, dv(:,1));
            dvd=map_sz2hts(h1, dv(:,2));
            dvd(dry) = nan;
            dud(dry) = nan;
            
            u=ob.xy.H*double(dud);
            v=ob.xy.H*double(dvd);
            
            u(badnodes,:) = nan;
            v(badnodes,:) = nan;
            
            
            % Construct vertical grid
            % Read timestep of elevations
            de=sz_readTimeStep(he,step);
            e=ob.xy.H*double(de);
            dp=ob.xy.H*gr.hgrid.depth;
            sz=sz_computeZlevels(dp,e,gr.vgrid);
            sz(badnodes,:) = nan;
            
            % calculate width of each layer
            wid = sz(:,2:end) - sz(:,1:end-1);
            
            % calculate vertically integrated s, u and average v
            av_u = (u(:,1:end-1) +  u(:,2:end) ).*wid./2;
            av_v = (v(:,1:end-1) +  v(:,2:end) ).*wid./2;
            av_var1 = (var1(:,1:end-1) +  var1(:,2:end) )./2;
            
            % calculate the total cross and along segment flow
            tot_u = 0.5 * (av_u(1:end-1,:) + av_u(2:end,:)).*trDist;
            tot_v = 0.5 * (av_v(1:end-1,:) + av_v(2:end,:)).*trDist;
            tot_var1 = 0.5 * (av_var1(1:end-1,:) + av_var1(2:end,:));
%            tot_var2 = 0.5 * (av_var2(1:end-1,:) + av_var2(2:end,:));
            
            area = trDist .* (wid(1:end-1,:) + wid(2:end,:))./2;
            
            tot_cross = tot_u .* x_comp + tot_v .* y_comp;
            tot_along = tot_u .* y_comp + tot_v .* -x_comp;
            
            tot_cross_var1 = tot_var1 .* tot_cross;
            
            in = tot_cross > 0;
            out = tot_cross < 0;
            indisp = tot_cross;
            outdisp = tot_cross;
            indisp(out) = indisp(out).*0;
            outdisp(in) = outdisp(in).*0;
            invar1 = tot_cross_var1;
            invar1(out) = invar1(out).*0;
            outvar1 = tot_cross_var1;
            outvar1(in) = outvar1(in).*0;
            
            cross(end+1) = nansum(tot_cross(:));
            along(end+1) = nansum(tot_along(:));
            
            var1inflow(end+1) = nansum(tot_cross_var1(in));
            var1outflow(end+1) = nansum(tot_cross_var1(out));
            
            inflow(end+1) = nansum(tot_cross(in));
            outflow(end+1) = nansum(tot_cross(out));
            
            % plots
            
            iw = 1024;
            ih = 400;
%            if rem(i,4)==1 && i < 90*7 && ii == 1
%                myfig = figure;
%                ph=pcolor(trLen(1:last,:), sz(1:last,2:end), tot_var1(1:last,:));hold on;
%                set(ph,'EdgeColor','none','FaceColor','interp');
%                axis([-inf inf ylims])
%                buf = sprintf('%s time step %d', runname, i);
%                title(buf);hold on;
%                cb = colorbar('EastOutside');
%                set(get(cb,'ylabel'),'string',var1name);
%                %caxis(cax); hold on;
%                
%                
%                fn = sprintf('%s-%s-%08d.png', runname, var1name, i * 960);
%                
%                set(gcf,'PaperUnits','inches','PaperPosition',[0 0 iw/100.0 ih/100.0])
%                print('-dpng', fn, '-r100');
%                close(myfig);
%                
%                % variable plot
%                myfig = figure;
%                ph=pcolor(trLen(1:last,:), sz(1:last,1:end-1), tot_cross_var1(1:last,:));hold on;
%                set(ph,'EdgeColor','none','FaceColor','interp');
%                axis([-inf inf ylims])
%                buf = sprintf('%s time step %d', runname, i);
%                title(buf);hold on;
%                cb = colorbar('EastOutside');
%                %caxis(cax); hold on;
%                set(get(cb,'ylabel'),'string',['cross flow (inflow positive) of ',var1name]);
%                
%                fn = sprintf('./%s-%s-%08d.png', runname, ['cross-',var1name], i * 960);
%                
%                set(gcf,'PaperUnits','inches','PaperPosition',[0 0 iw/100.0 ih/100.0])
%                print('-dpng', fn, '-r100');
%                close(myfig);
%                
%                % outflow var plot
%                myfig = figure;
%                ph=pcolor(trLen(1:last,:), sz(1:last,1:end-1), abs(outvar1(1:last,:)));hold on;
%                set(ph,'EdgeColor','none','FaceColor','interp');
%                axis([-inf inf ylims])
%                buf = sprintf('%s time step %d', runname, i);
%                title(buf);hold on;
%                cb = colorbar('EastOutside');
%                set(get(cb,'ylabel'),'string',[var1name,' outflow']);
%                
%                fn = sprintf('./%s-%s-%08d.png', runname, [var1name, '-outflow'], i * 960);
%                
%                set(gcf,'PaperUnits','inches','PaperPosition',[0 0 iw/100.0 ih/100.0])
%                print('-dpng', fn, '-r100');
%                close(myfig);
%                
%                % inflow var plot
%                myfig = figure;
%                ph=pcolor(trLen(1:last,:), sz(1:last,1:end-1), abs(invar1(1:last,:)));hold on;
%                set(ph,'EdgeColor','none','FaceColor','interp');
%                axis([-inf inf ylims])
%                buf = sprintf('%s time step %d', runname, i);
%                title(buf);hold on;
%                cb = colorbar('EastOutside');
%                set(get(cb,'ylabel'),'string',[var1name, 'inflow']);
%                
%                fn = sprintf('./%s-%s-%08d.png', runname, [var1name, '-inflow'], i * 960);
%                
%                set(gcf,'PaperUnits','inches','PaperPosition',[0 0 iw/100.0 ih/100.0])
%                print('-dpng', fn, '-r100');
%                close(myfig);
%                
%                % u plot
%                myfig = figure;
%                ph=pcolor(trLen(1:last,:), sz(1:last,1:end-1), tot_u(1:last,:)./area(1:last,:));hold on;
%                set(ph,'EdgeColor','none','FaceColor','interp');
%                axis([-inf inf ylims])
%                buf = sprintf('%s time step %d', runname, i);
%                title(buf);hold on;
%                cb = colorbar('EastOutside');
%                set(get(cb,'ylabel'),'string','E-W vel (m/s)');
%                caxis([-2 2]); hold on;
%                
%                
%                fn = sprintf('./%s-%s-%08d.png', runname, 'EW', i * 960);
%                
%                set(gcf,'PaperUnits','inches','PaperPosition',[0 0 iw/100.0 ih/100.0])
%                print('-dpng', fn, '-r100');
%                close(myfig);
%                
%                % v plot
%                myfig = figure;
%                ph=pcolor(trLen(1:last,:), sz(1:last,1:end-1), tot_v(1:last,:)./area(1:last,:));hold on;
%                set(ph,'EdgeColor','none','FaceColor','interp');
%                axis([-inf inf ylims])
%                buf = sprintf('%s time step %d', runname, i);
%                title(buf);hold on;
%                cb = colorbar('EastOutside');
%                caxis([-2 2]); hold on;
%                
%                
%                fn = sprintf('./%s-%s-%08d.png', runname, 'NS', i * 960);
%                
%                set(gcf,'PaperUnits','inches','PaperPosition',[0 0 iw/100.0 ih/100.0])
%                print('-dpng', fn, '-r100');
%                close(myfig);
%                
%                % cross plot
%                myfig = figure;
%                ph=pcolor(trLen(1:last,:), sz(1:last,1:end-1), tot_cross(1:last,:)./area(1:last,:));hold on;
%                set(ph,'EdgeColor','none','FaceColor','interp');
%                axis([-inf inf ylims])
%                buf = sprintf('%s time step %d', runname, i);
%                title(buf);hold on;
%                cb = colorbar('EastOutside');
%                set(get(cb,'ylabel'),'string','cross vel (inflow positive) (m/s)');
%                caxis([-2 2]); hold on;
%                
%                fn = sprintf('./%s-%s-%08d.png', runname, 'cross', i * 960);
%                
%                set(gcf,'PaperUnits','inches','PaperPosition',[0 0 iw/100.0 ih/100.0])
%                print('-dpng', fn, '-r100');
%                close(myfig);
%                
%                % along plot
%                myfig = figure;
%                ph=pcolor(trLen(1:last,:), sz(1:last,1:end-1), tot_along(1:last,:)./area(1:last,:));hold on;
%                set(ph,'EdgeColor','none','FaceColor','interp');
%                axis([-inf inf ylims])
%                buf = sprintf('%s time step %d', runname, i);
%                title(buf);hold on;
%                cb = colorbar('EastOutside');
%                set(get(cb,'ylabel'),'string','along boundary vel (m/s)');
%                caxis([-2 2]); hold on;
%                
%                fn = sprintf('./%s-%s-%08d.png', runname, 'along', i * 960);
%                
%                set(gcf,'PaperUnits','inches','PaperPosition',[0 0 iw/100.0 ih/100.0])
%                print('-dpng', fn, '-r100');
%                close(myfig);
%                
%                % outflow plot
%                myfig = figure;
%                ph=pcolor(trLen(1:last,:), sz(1:last,1:end-1), abs(outdisp(1:last,:)));hold on;
%                set(ph,'EdgeColor','none','FaceColor','interp');
%                axis([-inf inf ylims])
%                buf = sprintf('%s time step %d', runname, i);
%                title(buf);hold on;
%                cb = colorbar('EastOutside');
%                set(get(cb,'ylabel'),'string','outflow (m^3/s)');
%                
%                fn = sprintf('./%s-%s-%08d.png', runname, 'outflow', i * 960);
%                
%                set(gcf,'PaperUnits','inches','PaperPosition',[0 0 iw/100.0 ih/100.0])
%                print('-dpng', fn, '-r100');
%                close(myfig);
%                
%                % inflow plot
%                myfig = figure;
%                ph=pcolor(trLen(1:last,:), sz(1:last,1:end-1), abs(indisp(1:last,:)));hold on;
%                set(ph,'EdgeColor','none','FaceColor','interp');
%                axis([-inf inf ylims])
%                buf = sprintf('%s time step %d', runname, i);
%                title(buf);hold on;
%                cb = colorbar('EastOutside');
%                set(get(cb,'ylabel'),'string','inflow (m^3/s)');
%                
%                fn = sprintf('./%s-%s-%08d.png', runname, 'inflow', i * 900);
%                
%                set(gcf,'PaperUnits','inches','PaperPosition',[0 0 iw/100.0 ih/100.0])
%                print('-dpng', fn, '-r100');
%                close(myfig);
            end
        end
    end

    myfig = figure;
    subplot(2,1,1)
    plot(t,var1inflow);
    hold on;
    plot(t,-var1outflow,'r');
    title([var1name,' inflow (blue) and outflow (red)']);
    ylim([0 3.5e6])
    subplot(2,1,2)
    plot(t,inflow);
    hold on;
    plot(t,-outflow,'r');
    ylim([0 3.5e4])
    title('inflow (blue) and outflow (red)');
    fn = sprintf('./tseries%s-%s.png', runname,var1name);
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 iw/100.0 ih/100.0])
    print('-dpng', fn, '-r100');
   % close(myfig);
end
