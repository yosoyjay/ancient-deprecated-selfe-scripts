function modpar = ptparamread(parfile)

%  This function reads the parameters from param.in for use with other MATLAB scripts
%  to plot particleTracking results.  Modified to work with SELFE 3.x param.in files
%  
%  Input: param.in
%  Output: Parameters required for ptrack to be written in particle.bp
%

fid = fopen(parfile,'r');

fgets(fid);  % ipre !pre-processing flag (to output obs.out)
fgets(fid);  % ntracers
fgets(fid);  % nonhydro !write mode
fgets(fid);  % imm !0: without bed deformation; 1: with bed deformation (e.g., tsunami)

ics = fscanf(fid,'%i',1);  	% ics
slam0 = fscanf(fid,'%f',1); 	% center of projection, has different name in SELFE 3.x param.in
sfea0 = fscanf(fid,'%f',1); 	% "       "

fgets(fid);  % ibcc
fgets(fid);  % itransport
fgets(fid);  % nrampbc
fgets(fid);  % drampbc
fgets(fid);  % ihot2
fgets(fid);  % ihorcon
fgets(fid);  % ihdif
fgets(fid);  % ihdrag
fgets(fid);  % bfric
fgets(fid);  % ncor
fgets(fid);  % icst
fgets(fid);  % ic_elev
fgets(fid);  % ibcc_mean
fgets(fid);  % indvel
fgets(fid);  % shapiro
fgets(fid);  % rmaxvel
fgets(fid);  % velmin_btrack
fgets(fid);  % ihhat
fgets(fid);  % inunfl
fgets(fid);  % shapiro

h0 = fscanf(fid,'%f',1); fgets(fid);  % h0 

fgets(fid);  % thetai

rnday = fscanf(fid,'%i',1); fgets(fid);  % rnday

fgets(fid);  % nramp
fgets(fid);  % dramp

dt = fscanf(fid,'%i',1); fgets(fid);% dt

fgets(fid);  % slvr_output_spool
fgets(fid);  % mxitn
fgets(fid);  % tolerance
fgets(fid);  % nadv 
fgets(fid);  % dtb_max1
fgets(fid);  % dtb_max2
fgets(fid);  % inter_st
fgets(fid);  % inter_mom
fgets(fid);  % kr_co
fgets(fid);  % blend_internal
fgets(fid);  % blend_bnd
fgets(fid);  % iupwind_t
fgets(fid);  % nws
fgets(fid);  % wtiminc
fgets(fid);  % nrampwind
fgets(fid);  % drampwind
fgets(fid);  % iwindoff
fgets(fid);  % ihconsv
fgets(fid);  % isconsv
fgets(fid);  % itur
fgets(fid);  % turb_met
fgets(fid);  % turb_stab
fgets(fid);  % inu_st
fgets(fid);  % step_nu
fgets(fid);  % vnh1
fgets(fid);  % vnf1
fgets(fid);  % vnh2
fgets(fid);  % vnf2
fgets(fid);  % depth_zsigma
fgets(fid);  % s1_mxnbt
fgets(fid);  % s2_mxnbt
fgets(fid);  % iwrite

nspool = fscanf(fid,'%i',1); fgets(fid);% nspool 
ihfskip = fscanf(fid,'%i',1); fgets(fid);% ihfskip 

modpar.ics = ics;
modpar.h0 = h0;
modpar.rnday = rnday;
modpar.slam0 = slam0;
modpar.sfea0 = sfea0;
modpar.dt = dt;
modpar.nspool = nspool;
modpar.ihfskip = ihfskip;
