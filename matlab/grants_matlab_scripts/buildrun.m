%function out = buildrun(dummy,run_id,startday,endday,ncom_vers)
%out = buildrun(dummy,run_id,startday,endday,ncom_vers)
function [] = make_gr3(dummy)

% dummy is the path to the directory to make the gr3 files.

%ncom_vers = '/home/users/lawg/bin/readncomccs3.1_sirius2';
%app = '/home/users/lawg/bin/selfe/pelfe_3.1c_sirius_gotm';
%atmos = 3;% Use 1:NAM, 2:GFS, 3:NARR, 4:NNRP
estuary = true;
path(path, '/home/workspace/project/lopezj/grid_tests/grants_matlab_scripts/');



% Set up run directory
%[y1,m1,d1] = datevec(startday);
%[y2,m2,d2] = datevec(endday);
%jd0 = datenum(y1,1,1);
%wk1 = floor((startday-jd0)/7)+1;
%wk2 = floor((endday-jd0)/7);
%jd1 = jd0+(wk1-1)*7+1;
%jd2 = jd0+(wk2-1)*7+7;
%startday = datevec(jd1);
%endday = datevec(jd2);
%rundays = jd2-jd1+1;
%run_dir =  [num2str(y1) '-' num2str(wk1) '-' num2str(wk2) '-' run_id];
%eval(['!mkdir ' run_dir])
%eval(['!mkdir ' run_dir '/outputs'])
%eval(['!ln -sf ' dummy '/hgrid.gr3 ' run_dir '/'])
%eval(['!ln -sf ' dummy '/vgrid.in ' run_dir '/'])
%eval(['!ln -sf ' dummy '/param.in ' run_dir '/'])
eval(['cd ' dummy])
genhgridll;
if estuary
    genestgr3;
end
genalbgr3(.06,'constant');
gendraggr3([2e-3 5e-3],'depth');
gendiffgr3(1.0e-6,10,'constant');
gennudgegr3(5.2083e-4,'edge');
genintgr3;
genwtypegr3(7,'constant');
genxlscgr3(0.5,'constant');
%genbctides(jd1);
!LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH /home/users/lawg/bin/rotate_wind_spcs2ll -input hgrid.ll -output windrot_geo2proj.gr3 -ll2spcs
%
%eval(['!ln -sf ' app ' pelfe'])
%fg = hgrid2fg('hgrid.gr3');
%eval(['!/home/users/lawg/scripts/make_sflux_links.csh ' num2str(atmos) ' ' num2str(startday(1)) ' ' num2str(startday(2)) ' ' num2str(startday(3)) ' ' num2str(endday(1)) ' ' num2str(endday(2)) ' ' num2str(endday(3))])
%fid = fopen('run.qsub','w');
%fprintf(fid,'%s\n','#!/bin/bash');
%fprintf(fid,'%s\n','#$ -cwd');
%fprintf(fid,'%s\n','#$ -j y');
%fprintf(fid,'%s\n','#$ -S /bin/bash');
%%fprintf(fid,'%s\n','#$ -q eth.q'); %This line requests only eth0 nodes
%currdir = pwd;
%fprintf(fid,'%s\n',['/usr/mpi/gcc/openmpi-1.2.5/bin/mpirun --mca btl_tcp_if_include eth0 -np $NSLOTS ' currdir '/pelfe']);
%fclose(fid);
%
%if estuary
%    % Generate the temp.th file
%    wkstr = num2str(wk1,'%.2i');
%    eval(['!cat /home/workspace/ccalmr/elcirc/inputs/' num2str(y1) '/temp90_' num2str(y1) '_beaverfra.th_w' wkstr ' > tmp' wkstr '.dat'])
%    for i = (wk1+1):wk2
%        wkstr = num2str(i,'%.2i');
%        eval(['!cat tmp' num2str(i-1) '.dat /home/workspace/ccalmr/elcirc/inputs/' num2str(y1) '/temp90_' num2str(y1) '_beaverfra.th_w' wkstr ' > tmp' num2str(i) '.dat'])
%    end
%    eval(['!awk ''{print NR * 90,$2,$3}'' tmp' num2str(i) '.dat >temp.th']) 
%    !rm tmp*.dat
%
%    % Generate the flux.th file
%    wkstr = num2str(wk1,'%.2i');
%    eval(['!cat /home/workspace/ccalmr/elcirc/inputs/' num2str(y1) '/flux90_' num2str(y1) '_beaverfra.th_w' wkstr ' > tmp' wkstr '.dat'])
%    for i = (wk1+1):wk2
%        wkstr = num2str(i,'%.2i');
%        eval(['!cat tmp' num2str(i-1) '.dat /home/workspace/ccalmr/elcirc/inputs/' num2str(y1) '/flux90_' num2str(y1) '_beaverfra.th_w' wkstr ' > tmp' num2str(i) '.dat'])
%    end
%    eval(['!awk ''{print NR * 90,$2,$3}'' tmp' num2str(i) '.dat >flux.th']) 
%    !rm tmp*.dat
%end
%
%% Update the vgrid.in file. The example vgrid file must have the same
%% number of Z-levels listed in the first line. All other parameters can be
%% changed in the original file, and the new file will have updated values.
%% The Z-levels will be spaced smoothly with the last few S-levels.
%fid1 = fopen('vgrid.in');
%fid2 = fopen('vgrid.in.new','w');
%out = fscanf(fid1,'%i %i %f\n',3);zlev=out(2,1);slev=out(1,1)-zlev+1;sztrans=out(3,1);
%fgets(fid1); fgets(fid1);
%out = fscanf(fid1,'%i %f\n',2);z=out(2,1);
%out = fgets(fid1);
%while ~any(strfind(out,'S'))
%    out = fgets(fid1);
%end
%out = fscanf(fid1,'%f %f %f\n',3);hc=out(1,1);thetab=out(2,1);thetaf=out(3,1);
%% Now calculate all the levels
%sigma = (-1:1/(slev-1):0);
%if zlev==1
%    zdepth = z;
%else
%    %if thetab>0
%    %    thetab = 70/sztrans; %Place deep constriction at pycnocline
%    %end
%    zdat = sigmacalc(hc,thetab,thetaf,slev,sztrans,z);
%    last = zdat(slev-1:slev+1,end); rat = diff(last(2:3))/diff(last(1:2));
%    zdepth = -sztrans;
%    cnt = 1;
%    zdiff = last(1)-last(3);
%    while (zdepth(cnt)-z)/zdiff > zlev-cnt
%        zdepth(cnt+1) = zdepth(cnt)-rat*zdiff;
%        zdiff = zdepth(cnt)-zdepth(cnt+1);
%        cnt = cnt+1;
%    end
%    zdiff = (zdepth(cnt)-z)/(zlev-cnt);
%    zdepth(cnt:zlev) = zdepth(cnt):-zdiff:z;
%end
%zdepth = fliplr(zdepth);
%fprintf(fid2,'%i %i %f\n',[zlev+slev-1 zlev sztrans]);
%fprintf(fid2,'%s\n','Z levels');
%for i = 1:zlev
%    fprintf(fid2,'%i %f\n',[i zdepth(i)]);
%end
%fprintf(fid2,'%s\n','S levels');
%fprintf(fid2,'%f %f %f\n',[hc thetab thetaf]);
%for i = 1:slev
%    fprintf(fid2,'%i %f\n',[i sigma(i)]);
%end
%fclose(fid1);
%fclose(fid2);
%!rm vgrid.in
%!mv vgrid.in.new vgrid.in
%
%% Determine which boundaries are forced by ocean
%obs = [];
%for i = 1:length(fg.ob)
%    if strcmp(fg.ob(i).type,'Ocean')
%        obs = [obs i];
%    end
%end
%
%% Generate the hotstart.in file
%fid1 = fopen('date.in','w');
%fid2 = fopen([dummy '/date.in'],'r');
%out = fscanf(fid2,'%i %i %i\n',3);
%fprintf(fid1,'%s\n',num2str([1 length(obs) obs]));
%fprintf(fid1,'%i %i %i %i\n',[y1 m1 d1 0]);
%fprintf(fid1,'%i\n',0);
%fprintf(fid1,'%s\n','Selfe');
%fprintf(fid1,'%f %f\n',[10 33]);
%fprintf(fid1,'%i\n',0);
%fclose(fid1);
%fclose(fid2);
%eval(['!' ncom_vers])
%
%% Generate the temp_nu.in file
%fid1 = fopen('date.in','w');
%fid2 = fopen([dummy '/date.in'],'r');
%out = fscanf(fid2,'%i %i %i\n',3);
%fprintf(fid1,'%s\n',num2str([1 length(obs) obs]));
%fprintf(fid1,'%i %i %i %i\n',[y1 m1 d1 rundays]);
%fprintf(fid1,'%i\n',0);
%fprintf(fid1,'%s\n','Selfe');
%fprintf(fid1,'%f %f\n',[10 33]);
%fprintf(fid1,'%i\n',0);
%fclose(fid1);
%fclose(fid2);
%eval(['!' ncom_vers])
%
% Generate the uv_bcc.th file
%fid1 = fopen('date.in','w');
%fid2 = fopen([dummy '/date.in'],'r');
%out = fscanf(fid2,'%i %i %i\n',3);
%fprintf(fid1,'%s\n',num2str([2 length(obs) obs]));
%fprintf(fid1,'%i %i %i %i\n',[y1 m1 d1 rundays]);
%fprintf(fid1,'%i\n',0);
%fprintf(fid1,'%s\n','Selfe');
%fprintf(fid1,'%f %f\n',[10 33]);
%fprintf(fid1,'%i\n',0);
%fclose(fid1);
%fclose(fid2);
%eval(['!' ncom_vers])
%
%% No uv3D.th file being generated 
%
%% Update rnday variable in param.in
%fid1 = fopen('param.in');
%fid2 = fopen('param.in.new','w');
%out = fgets(fid1);
%while out~=-1
%    if strcmp(out,'rnday')
%        fprintf(fid2,'%s',['  rnday = ' num2str(rundays)]);
%    else
%        fprintf(fid2,'%s',out);
%    end
%    out = fgets(fid1);
%end
%fclose(fid1);
%fclose(fid2);
%!rm param.in
%!mv param.in.new param.in
%
