function modpar = ptparamread(parfile)

fid = fopen(parfile,'r');
fgets(fid);  % version number
fgets(fid);  % start time
fgets(fid);  % ipre !pre-processing flag (to output obs.out)
fgets(fid);  % nscreen
fgets(fid);  % iwrite !write mode
fgets(fid);  % imm !0: without bed deformation; 1: with bed deformation (e.g., tsunami)
fgets(fid);  % ihot
ics = fscanf(fid,'%i',1); fgets(fid);  % ics

out = fscanf(fid,'%f %f',2); fgets(fid);% slam0,sfea0
slam0 = out(1); sfea0 = out(2);

fgets(fid);  % thetai
fgets(fid);  % ibc,ibtp
fgets(fid);  % nrampbc,drampbc 
rnday = fscanf(fid,'%i',1); fgets(fid);  % rnday
fgets(fid);  % nramp,dramp
dt = fscanf(fid,'%i',1); fgets(fid);% dt
fgets(fid);  % !nsubfl !flag
fgets(fid);  % nadv !flag: 1-Euler; 2: R-K
fgets(fid);  % dtb_max1,dtb_max2 !min. or max. sub-step for Euler and R-K if nadv=0; otherwise only dtb_max1 is used
h0 = fscanf(fid,'%f',1); fgets(fid);  % h0 
fgets(fid);  % nchi
fgets(fid);  % ncor
fgets(fid);  % nws,wtiminc
fgets(fid);  % coldstart
fgets(fid);  % windoff
fgets(fid);  % ihconsv,isconsv
fgets(fid);  % itur
fgets(fid);  % mid,stab
fgets(fid);  %icst
fgets(fid);  % ntip,tip_dp !cut-off depth for applying tidal potential

nbfr =fscanf(fid,'%i',1); fgets(fid); % nbfr

for i=1:nbfr
  fgets(fid);  % !tag
  fgets(fid);  % amig(i),ff(i),face(i) !freq., nodal factor and earth equil.
end

nope=fscanf(fid,'%i\n',1); fgets(fid); % nope

out = ones(1,nope);
ntmp=out; iettype=out; ifltype=out; itetype=out; isatype=out;
for k=1:nope
  out = fscanf(fid,'%i %i %i %i %i ',5); fgets(fid); % ntmp,iettype(k),ifltype(k),itetype(k),isatype(k),ifitype(k)
  ntmp(k)=out(1); iettype(k)=out(2); ifltype(k)=out(3); itetype(k)=out(4); isatype(k)=out(5);
  if iettype(k)==2
    fgets(fid);  % eth(k,1)
  else if iettype(k)==3
    for i=1:nbfr
      fgets(fid);  %  !freq. name
      for j=1:ntmp(k)
        fgets(fid);  % emo(k,j,i),efa(k,j,i) !amp. and phase
      end
    end
      end
  end

  if ifltype(k)~=0
  if ifltype(k)==3
    for i=1:nbfr
      fgets(fid);  %
      fgets(fid);  % vmo(k,i),vfa(k,i) !uniform amp. and phase along each segment
    end
  end
  
  if itetype(k)==2
    fgets(fid);  % tth(k,1,1)
  else
      if itetype(k)==-1
          fgets(fid);  % tobc !nudging factor
      end
  end

  if isatype(k)==2
    fgets(fid);  % sth(k,1,1)
  else
      if isatype(k)==-1
          fgets(fid);  % sobc !nudging factor
      end
  end
  end

%  if ifitype(k)==2
%    out = fgets(fid);  % fth(k,1,1)
%  else if ifitype(k)==-1
%    out = fgets(fid);  % fobc !nudging factor
%      end
%  end
end

fgets(fid); 
out = fscanf(fid,'%i %i\n',2); % nspool,ihfskip !output and file spools
nspool = out(1);
ihfskip = out(2);

modpar.ics = ics;
modpar.h0 = h0;
modpar.rnday = rnday;
modpar.slam0 = slam0;
modpar.sfea0 = sfea0;
modpar.dt = dt;
modpar.nspool = nspool;
modpar.ihfskip = ihfskip;