function status = genbctides(jd)
% 
[yr1,m1,d1] = datevec(jd);
fg1 = hgrid2fg('hgrid.gr3');
fg2 = hgrid2fg('/home/users/yinglong/AMB24Stuff/ElcircScripts/TIDES2/edpac2xy.gr3');
fid = fopen('/home/users/yinglong/AMB24Stuff/ElcircScripts/TIDES2/teanl.tct');
fnames = {'Z0' 'O1' 'K1' 'Q1' 'P1' 'K2' 'N2' 'M2' 'S2'};
freqs = length(fnames);
fgets(fid);
out = fscanf(fid,'%i %f %i\n',3); nodes = out(3);
tides.x =  fg2.x;
tides.y =  fg2.y;
ob_ind = [];
totb = length(fg1.ob);
for i = 1:totb
    if strcmp(fg1.ob(i).type,'Ocean')
        ob_ind = [ob_ind i];
    end
end
obs = length(ob_ind);
el_flag = zeros(1,obs);
vel_flag = zeros(1,obs);
t_flag = zeros(1,obs);
s_flag = zeros(1,obs);

for i = 1:9 % Only read in the first 9 frequencies
    tides.data(i).forcing = fscanf(fid,'%f\n',1); fgets(fid);
    fgets(fid);% Frequency names aren't in the teanl.tct file for some reason...
    data = fscanf(fid,'%i %f %f\n',[3,nodes]);
    tides.data(i).amplitude = data(2,:)';
    tides.data(i).phase =  data(3,:)';
    tides.data(i).name =  fnames{i};
end

% Convert x/y to same format as tidal model grid
[fg1.x,fg1.y] = convm2ll(fg1.x,fg1.y);
[fg1.x,fg1.y] = convll2m(fg1.x,fg1.y,'nos8');

% Generate list of boundary nodes to interpolate
for b = 1:totb
    if strcmp(fg1.ob(b).type,'Ocean')
        bnodes(b).ind = fg1.ob(b).ind';
        el_flag(b) = 3;% Forced by tides
        vel_flag(b) = 0;% Normal vel. not specified
        t_flag(b) = 0;% Temperature not specified
        s_flag(b) = 0;% Salinity not specified
    else
        bnodes(b).ind = [];
        el_flag(b) = 0;% Not forced by tides
        vel_flag(b) = 1;% Normal vel. not specified
        t_flag(b) = 1;% Temperature not specified
        s_flag(b) = 2;% Constant salinity (sthconst=0.0)
    end
end

% Generate nearest-nodes list
for b = 1:totb
    if strcmp(fg1.ob(b).type,'Ocean')
        for j = 1:length(fg1.ob(b).ind)
            dist = sqrt((tides.x-fg1.x(fg1.ob(b).ind(j))).^2 + (tides.y-fg1.y(fg1.ob(b).ind(j))).^2);
            [dist,ind] = sort(dist);
            bnodes(b).wind(j,:) = ind(1:3);
        end
    end
end

% Use interpolation method from intel_deg.f
for b = 1:totb
    if strcmp(fg1.ob(b).type,'Ocean')
        x = fg1.x(bnodes(b).ind);
        y = fg1.x(bnodes(b).ind);
        n1 = bnodes(b).wind(:,1);
        n2 = bnodes(b).wind(:,2);
        n3 = bnodes(b).wind(:,3);
        x1 = fg2.x(n1);
        x2 = fg2.x(n2);
        x3 = fg2.x(n3);
        y1 = fg2.y(n1);
        y2 = fg2.y(n2);
        y3 = fg2.y(n3);
        a1 = x2.*y3 - x3.*y2;
        b1 = y2 - y3;
        c1 = x3 - x2;
        a2 = x3.*y1 - x1.*y3;
        b2 = y3 - y1;
        c2 = x1 - x3;
        a3 = x1.*y2 - x2.*y1;
        b3 = y1 - y2;
        c3 = x2 - x1;
        a = a1 + a2 + a3;
        rl1 = (a1+b1.*x+c1.*y)./a;
        rl2 = (a2+b2.*x+c2.*y)./a;
        rl3 = (a3+b3.*x+c3.*y)./a;
        clear elmod phase
        for freq = 1:freqs
            elmod(:,1) = tides.data(freq).amplitude(n1);
            elmod(:,2) = tides.data(freq).amplitude(n2);
            elmod(:,3) = tides.data(freq).amplitude(n3);
            phase(:,1) = tides.data(freq).phase(n1);
            phase(:,2) = tides.data(freq).phase(n2);
            phase(:,3) = tides.data(freq).phase(n3);
            atmp = rl1.*elmod(:,1).*cos(phase(:,1)) + rl2.*elmod(:,2).*cos(phase(:,2)) + rl3.*elmod(:,3).*cos(phase(:,3));
            btmp = rl1.*elmod(:,1).*sin(phase(:,1)) + rl2.*elmod(:,2).*sin(phase(:,2)) + rl3.*elmod(:,3).*sin(phase(:,3)); 
            phn = atan2(btmp,atmp);
            an = atmp./cos(phn);
            phn(phn<0.0) = phn(phn<0) + 2.0.*pi;
            phn = phn.*180./pi;
            bnodes(b).freq(freq).amplitude = an;
            bnodes(b).freq(freq).phase = phn;
        end
    end
end

% Generate the Z0 values for bctides.in
% (Invalid entries are not currently being corrected)
fid1 = fopen('date.in','w');
fprintf(fid1,'%i %i %i %i\n',[yr1 m1 d1 0]);
fprintf(fid1,'%i\n',0);
fprintf(fid1,'%s\n',num2str([obs ob_ind]));
fclose(fid1);
!/home/users/lawg/bin/readssh2b_sirius
fid = fopen('Z0.out');
for b = 1:totb
    if strcmp(fg1.ob(b).type,'Ocean')
        fgets(fid);
        out = zeros(length(bnodes(b).ind),1);
        for i = 1:length(bnodes(b).ind)
            out(i) = fscanf(fid,'%e',1); fgets(fid);
        end
        bnodes(b).freq(1).amplitude = out;
        bnodes(b).freq(1).phase = zeros(size(out));
    else
        bnodes(b).freq(1).amplitude = 0;
        bnodes(b).freq(1).phase = 0;
    end
end
!rm Z0.out

% Generate nodal factors and Earth equilibrium arg's 
fid = fopen('tide.in','w');
year1 = num2str(yr1); year2 = year1(3:4); year1 = year1(1:2);
fprintf(fid,'%s\n',['  08' num2str(d1,'%.2i') num2str(m1,'%.2i') year2 year1]);
fprintf(fid,'%s\n','     8615  SVI tidal model      GMT 461512400');
fclose(fid);
!cp /home/users/lawg/matlab/my_collection/tides/fort.8 .
!/home/users/lawg/bin/./tide1_e < tide.in > tide.out
!rm fort.8 
for i = 1:freqs
    fid = fopen('tide.out');
    out = '0000000000';
    while ~strcmp(strtrim(out(1:8)),fnames{i})&out(1)~=-1
        out = fgets(fid);
    end
    if out(1)==-1
        error(['Couldn''t find the ' fnames{i} ' frequency in the tide.out file'])
    end
    tides.data(i).nf = str2num(out(8:end));
end
!rm tide.out
!rm tide.in

% Now write the bctides.in file
% Set flags
fid = fopen('bctides.in','w');
fprintf(fid,'%s\n',datestr(jd,0));
fprintf(fid,'%s\n','0 40 ntip');% [# of constituents in earth tidal potential   depth above which not to calculate tides]
fprintf(fid,'%s\n',[num2str(freqs,'%i') ' nbfr']);
for f = 1:freqs
    fprintf(fid,'%s\n',fnames{f});
    fprintf(fid,'%e %f %f\n',[tides.data(f).forcing tides.data(f).nf]);
end
fprintf(fid,'%s\n',[num2str(length(fg1.ob)) ' nope']);
for b = 1:totb
    fprintf(fid,'%i %i %i %i %i\n',[length(fg1.ob(b).ind) el_flag(b) vel_flag(b) t_flag(b) s_flag(b)]);
    if strcmp(fg1.ob(b).type,'Ocean')
        for f = 1:freqs
            fprintf(fid,'%s\n',fnames{f});
            fprintf(fid,'%e %f\n',[bnodes(b).freq(f).amplitude bnodes(b).freq(f).phase]');
        end
    else
        fprintf(fid,'%f\n',0);
    end
end
fclose(fid);