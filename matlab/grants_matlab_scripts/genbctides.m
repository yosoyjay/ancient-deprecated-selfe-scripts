function status = genbctides(jd)
% 
method = 'standard';% Use 'notide' or 'standard'

[yr1,m1,d1] = datevec(jd);
fg1 = hgrid2fg('hgrid.gr3');
fnames = {'Z0' 'O1' 'K1' 'Q1' 'P1' 'K2' 'N2' 'M2' 'S2'};
freqs = length(fnames);
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

% Build nos8 version of hgrid .gr3
fid8 = fopen('hgrid.nos8','w');
fprintf(fid8,'%s\n','Nos8 version of hgrid for intel_deg input');
fprintf(fid8,'%i %i\n',[length(fg1.e) length(fg1.x)]);
[lon,lat] = convm2ll(fg1.x,fg1.y);
[fg1.x,fg1.y] = convll2m(lon,lat,'nos8');
block = [1:length(fg1.x);fg1.x';fg1.y';fg1.z'];
fprintf(fid8,'%i %f %f %f\n',block);
block = [1:length(fg1.e);ones(1,length(fg1.e)).*3;fg1.e'];
fprintf(fid8,'%i %i %i %i %i\n',block);
fid9 = fopen('hgrid.gr3');
for i = 1:length(fg1.x)+length(fg1.e)+2
    fgets(fid9);
end
out = fgets(fid9);
while out~=-1
    fprintf(fid8,'%s\n',out);
    out = fgets(fid9);
end
fclose(fid8);
fclose(fid9); 

% Build station files for intel_deg to chomp on...
for b = 1:totb
    if strcmp(fg1.ob(b).type,'Ocean')
        bnodes(b).ind = fg1.ob(b).ind';
        el_flag(b) = 3;% Forced by tides
        vel_flag(b) = 0;% Normal vel. not specified
        t_flag(b) = 0;% Temperature not specified
        s_flag(b) = 0;% Salinity not specified
        fid(b).id = fopen(['boundary' num2str(b) '.sta'],'w');
        fprintf(fid(b).id,'%s\n','Station file for input to intel_deg');
        fprintf(fid(b).id,'%i\n',length(fg1.ob(b).ind));
        fprintf(fid(b).id,'%i %f %f\n',[(1:length(fg1.ob(b).ind));fg1.x(fg1.ob(b).ind)';fg1.y(fg1.ob(b).ind)']);
        fclose(fid(b).id);
    else
        bnodes(b).ind = [];
        el_flag(b) = 0;% Not forced by tides
        vel_flag(b) = 1;% Normal vel. not specified
        t_flag(b) = 1;% Temperature not specified
        s_flag(b) = 2;% Constant salinity (sthconst=0.0)
    end
end
% Do the calkoolatin'
for i = 1:obs
    if strcmp(fg1.ob(i).type,'Ocean')
        eval(['!mv boundary' num2str(i) '.sta intel_deg.sta'])
        !/home/users/lawg/bin/intel_deg
        fid = fopen('intel_deg.out');
        fgets(fid);
        bns = fscanf(fid,'%i\n',1);
        fgets(fid);
        for j = 1:9
            tides.data(j).init = fscanf(fid,'%f\n',1);
            fgets(fid);
            out = fscanf(fid,'%f %f\n',[2,bns]);
            bnodes(i).freq(j).amplitude = out(1,:)';
            bnodes(i).freq(j).phase = out(2,:)';
            bnodes(i).freq(j).name = fnames{j};
        end
    end
end
!rm intel_deg.sta

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
    switch true
        case strcmp(method,'standard')
            fprintf(fid,'%e %f %f\n',[tides.data(f).init tides.data(f).nf]);
        case strcmp(method,'notide')
            fprintf(fid,'%e %f %f\n',[tides.data(f).init.*0 tides.data(f).nf.*0]);
    end
end
fprintf(fid,'%s\n',[num2str(length(fg1.ob)) ' nope']);
for b = 1:totb
    fprintf(fid,'%i %i %i %i %i\n',[length(fg1.ob(b).ind) el_flag(b) vel_flag(b) t_flag(b) s_flag(b)]);
    if strcmp(fg1.ob(b).type,'Ocean')
        for f = 1:freqs
            fprintf(fid,'%s\n',fnames{f});
            switch true
                case strcmp(method,'standard')
                    fprintf(fid,'%e %f\n',[bnodes(b).freq(f).amplitude bnodes(b).freq(f).phase]');
                case strcmp(method,'notide')
                    fprintf(fid,'%e %f\n',[bnodes(b).freq(f).amplitude.*0 bnodes(b).freq(f).phase.*0]');
            end
        end
    else
        fprintf(fid,'%f\n',0);
    end
end
fclose(fid);
!rm fort.* date.in intel_deg.out hgrid.nos8