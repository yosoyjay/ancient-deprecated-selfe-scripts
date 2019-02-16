function status = geninitcond(method,srange,trange)
% METHOD is either 'vertical' or 'horizontal'
% SRANGE is either one value (for vertically-constant salinities)
%       or two values in the form [minval maxval]
% TRANGE is either one value (for vertically-constant temperatures)
%       or two values in the form [minval maxval]

vmeth = 'linear';% 'linear' or 'sigmoid' Need to make this a feature of the function
if srange(1)>srange(2); srange = fliplr(srange); end
if trange(1)>trange(2); trange = fliplr(trange); end

switch method
    case 'vertical'
        fid = fopen('vgrid.in');
        out = fscanf(fid,'%i %i %f',3);levs=out(1);slevs=out(1)-out(2)+1;zlevs=out(2);sztrans=out(3);fgets(fid);
        out = fgets(fid);
        for i = 1:zlevs
            out = fscanf(fid,'%i %f\n',2);
            zlev(i) = -out(2);
        end
        fg = hgrid2fg('hgrid.gr3');
        fgets(fid);
        out = fscanf(fid,'%f %f %f\n',3);contrans=out(1);pos=out(2);constr=out(3);
        for i = 1:slevs
            out = fscanf(fid,'%i %f\n',2);
            slev(i) = out(2);
        end
        zdepth = sigmacalc(contrans,pos,constr,slevs,sztrans,zlev(1));
        zdepth = [-zdepth(end,1:end-2)';flipud(zlev')];zdepth(1) = 0;
        zmax = zdepth(end);
        
        if length(trange)==1
            t = ones(size(zdepth)).* trange;
        else
            switch vmeth
                case 'linear'
                    mt = -diff(trange)/zdepth(end-4);
                    t = [(zdepth(1:end-3).*mt+trange(2));trange(1).*ones(3,1)];
                case 'sigmoid'
                    if pos>0&pos<1
                        zcline = pos.*sztrans;
                    else
                        if sztrans<100
                            zcline = -sztrans;
                        else
                            zcline = sztrans/-2;
                        end
                    end
                    b = .2;
                    c = zcline./5;
                    if length(trange)==2
                        a = (trange(2)-trange(1))./pi;
                        d = mean(trange);
                        t = a.*atan((b.*zdepth)+c)+d;
                    end
            end
        end
        if length(srange)==1
            s = ones(size(zdepth)).* srange;
        else
            switch vmeth
                case 'linear'
                    ms = diff(srange)/zdepth(end-4);
                    s = [(zdepth(1:end-3).*ms + srange(1));srange(2).*ones(3,1)];
                case 'sigmoid'
                    if pos>0&pos<1
                        zcline = pos.*sztrans;
                    else
                        if sztrans<100
                            zcline = -sztrans;
                        else
                            zcline = sztrans/-2;
                        end
                    end
                    b = .2;
                    c = zcline./5;
                    if length(srange)==2
                        a = (srange(2)-srange(1))./pi;
                        d = mean(srange);
                        s = a.*atan((b.*zdepth)+c)+d;
                    end
            end
        end
        fid = fopen('ts.ic','w');
        fprintf(fid,'%s\n',[num2str(levs) ' ! Total number of levels']);
        block = [(1:levs);fliplr(-zdepth');fliplr(t');fliplr(s')];
        fprintf(fid,'%i %f %f %f\n',block);
        fclose(fid);
        status = exist('ts.ic','file');
    case 'horizontal'
    otherwise
        error('Didn''t recognize that method...');
end