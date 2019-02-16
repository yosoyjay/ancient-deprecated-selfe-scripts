function particles = layparticles(gridfile,modpar,partnum,relpos,reltime,clustsz)
% particles = LAYPARTICLES(gridfile,modpar,partnum,relpos,reltime);
% 
% gridfile:  Name of domain grid file
% modpar: Important parameters used in forcast model (use PTPARAMREAD)
% partnum:  Number of particles to release 
% relpos:  coordinates defining the points from which particles will
% be released
% relform:  Either 'Random', or 'Regular' distributions
% reltime:  The time of release (seconds from beginning of run)
%

reldepth = -0.5;
modpar.rnday = 1;% Force only 1 day worth of drifter tracking
relform = 'Regular';


if nargin<4
    reltime = 1;
    %relform = [-124 46.16 ;-124 46.29;-124.14 46.29;-124.14 46.16];
    relpos = [2.6390e+05 3.7661e+05;3.2384e+05 3.7798e+05];
end

switch 1
    case size(relpos,1)==1
        xpos = relpos(1,1);
        ypos = relpos(1,2);
        data = [1 reltime(1) xpos ypos reldepth];
    case size(relpos,1)==2
%    [trans(:,1),trans(:,2)] = convll2m(relpos(:,1),relpos(:,2));
        trans = relpos;
        x = trans(2,1) - trans(1,1);
        m = (trans(2,2)-trans(1,2))/(trans(2,1)-trans(1,1));
        b = trans(2,2)-trans(2,1)*m;
        if strcmp(relform,'Random')%random
            xpos = rand(partnum,1).*x + trans(1,1);
        else
            if strcmp(relform,'Regular') % Evenly spaced
                xpos = (trans(1,1):(x/(partnum-1)):trans(2,1))';
            end
        end
        ypos = (m.*xpos + b);
        data = [(1:partnum)' ones(partnum,1).*reltime(1) xpos ypos ones(partnum,1).*reldepth];
    case size(relpos,1)>=3
        % Create release clusters centered on the relpos positions
        clustnum = partnum/size(relpos,1);
        dx = sqrt((pi*clustsz.big^2)/clustsz.big);
        drops = size(relpos,1);
        brad = clustsz.big;% Includes all particles released within clustsz.big meters
        bpnts = ceil(brad/10);
        brho=ones(1,bpnts)*brad;
        btheta=linspace(0,2*pi,bpnts);
        [bcirx,bciry] = pol2cart(btheta,brho);
        xtot = [];
        ytot = [];
        for i = 1:size(relpos,1)
            % Estimate relpos using a square
            xmin = relpos(i,1)-clustsz.big;
            xmax = relpos(i,1)+clustsz.big;
            xlen = xmax-xmin;
            ymin = relpos(i,2)-clustsz.big;
            ymax = relpos(i,2)+clustsz.big;
            ylen = ymax-ymin;
            if strcmp(relform,'Regular')
                [xpos,ypos] = meshgrid((xmin:dx:xmax),(ymin:dx:ymax));
                x = xpos(:);
                y = ypos(:);
            else
                x = rand(clustnum*1.3,1).*xlen + xmin;
                y = rand(clustnum*1.3,1).*ylen + ymin;
            end
            
            % Find particles lying outside the domain
            
            bpolx=bcirx+relpos(i,1);
            bpoly=bciry+relpos(i,2);
            
            inside = inpolygon(x,y,bpolx,bpoly);
            x(~inside) = [];
            y(~inside) =[];
            % Find particles lying on islands
            isles = sum(isnan(gridfile.bndpth))-2;
            isle = find(isnan(gridfile.bndpth)); isle(1) = [];
            for j = 1:isles
                islearray = gridfile.bndpth(isle(j)+1:isle(j+1)-1);
                dry = inpolygon(x,y,gridfile.x(islearray),gridfile.y(islearray));
                x(dry) = [];
                y(dry) =[];
            end
            xtot = [xtot;x];
            ytot = [ytot;y];
        end
        x = xtot;
        y = ytot;
        num = length(x);
        data = [(1:num)' ones(num,1).*reltime(1) x y ones(num,1).*reldepth];
        
end
   
% Now set up additional releases
nrel = length(reltime);
if nrel>1
    block = zeros(num*(nrel-1),5);
    for i = 1:nrel-1
        block((i-1)*num+1:i*num,:) = data;
        block((i-1)*num+1:i*num,2) = ones(num,1).*reltime(i+1);
    end
    partnum = num*nrel;
    data = [data;block];
    data(:,1) = (1:partnum);
end

particles = [xpos ypos];

fid = fopen('particle.bp','w');
fprintf(fid,'%s\n','Particle set for drifter tracking forecast');
fprintf(fid,'%i\n',1);% 1 = Write to screen
fprintf(fid,'%i\n',1);% ibf=1 for forward tracking
fprintf(fid,'%i\n',1);% 1 = Constrain particle depth to initial condition
fprintf(fid,'%i %5.2f %4.2f\n',[1 modpar.slam0 modpar.sfea0]);
fprintf(fid,'%2.1f %i %i %i %i %i\n',[modpar.h0 modpar.rnday modpar.dt modpar.nspool modpar.ihfskip 10]);
fprintf(fid,'%i\n',partnum);
fprintf(fid,'%5i %6i %8.2f %8.2f %4.2f\n',data');
fclose(fid);
