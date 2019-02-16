function out = genvgrid(slevs,zlevs,sztrans,constr,pos,contrans)
%GENVGRID Utility for generating vgrid.in files for SELFE runs.
%   OUT = GENVGRID(slevs,zlevs,sztrans,constriction,depth,contrans),
%   Using the required inputs, GENVGRID creates a vgrid.in file to be used
%   in a SELFE circulation model run. The function returns a value of 1 if
%   it successfully creates the file, and a 0 if not.
%  
%   SLEVS is the number of sigma-levels desired
%   ZLEVS is the number of Z-levels desired
%   SZTRANS is the depth at which vertical levels should transition from
%       sigma to z-layer forms 
%   CONSTR is the factor by which sigma-levels constrict at the surface and
%       subsurface positions
%   POS defines the vertical location of the subsurface constriction. A
%       value of 1.0 will place the subsurface constriction at the depth of
%       SZTRANS (not recommended if any Z-levels are being used). 0.5 
%       places the constriction midway between SZTRANS and the surface.
%       0.0 places the constriction at the surface (because a constriction
%       already exists at the surface, this option essentially eliminates
%       the subsurface constriction).
%   CONTRANS is the isobath at which the vertical distribution of sigma
%       levels transitions from an even distribution, to one defined by the
%       constriction factors CONSTR and DEPTH.
%
%   EXAMPLE:
%      out = genvgrid(30,20,100,5,0,30);
%

%-------------------------------
%   Additional details:
%
%   C. Grant Law 04-04-2011

if sztrans<0.0; sztrans = -sztrans; end
if contrans<0.0; contrans = -contrans; end
if pos==1.0 & zlevs>0 & slevs>0
    error(['When performing a run with hybrid SZ-levels, the value for POS should not be 1.0. Instead, either reposition the SZ-transition to a deeper depth and set the POS value to place the subsurface constriction where you want it (i.e., POS = old_sztrans/new_sztrans), or set POS to 1.0 to eliminate the subsurface constriction entirely.'])
end
if ~exist('hgrid.gr3','file')
    error('A valid hgrid.gr3 file is needed in the current directory')
else
    fg = hgrid2fg('hgrid.gr3');
    depth = -(min(fg.z)) + 1;
end
if zlevs==0;
    zlevs = 1;
end
if slevs==0; slevs = 10; sztrans = 10; constr = 0.003; pos = 0.0; contrans = 9.9; end% This optimizes for z-level only runs

sigma = (-1:1/(slevs-1):0);
if zlevs==1
    zdepth = depth;
else
    if pos == 0
        zdat = sigmacalc(contrans,pos,constr,slevs+zlevs,depth,depth);
        zdepth = -zdat(end,slevs:end-1);
        zlevs = length(zdepth);
        sztrans = zdepth(1);
    else
        zdat = sigmacalc(contrans,pos,constr,slevs,sztrans,depth);
        last = zdat(end,end-3:end-1);
        rat = max([1.1 diff(last(2:3))/diff(last(1:2))]);
        zdepth = sztrans;
        cnt = 1;
        zdiff = last(3)-last(2);
        while (((zdepth(cnt)-depth)/zdiff > zlevs-cnt) & cnt<=zlevs)
            zdepth(cnt+1) = zdepth(cnt)-rat*zdiff;
            zdiff = zdepth(cnt)-zdepth(cnt+1);
            cnt = cnt+1;
        end
        if cnt>zlevs
            short = (depth-zdepth(end))./(zdepth(end)-zdepth(end-1));
            error(['Need about ' num2str(ceil(short)) ' more Z-levels to maintain a constant rate of increase in layer thicknesses'])
            zdepth = zdepth(1:zlevs);
        else
            zdiff = (zdepth(cnt)-depth)/(zlevs-cnt);
            zdepth(cnt:zlevs) = zdepth(cnt):-zdiff:depth;
        end
    end
end
zdepth = fliplr(zdepth);
fid = fopen('vgrid.in','w');
fprintf(fid,'%i %i %f\n',[zlevs+slevs-1 zlevs sztrans]);
fprintf(fid,'%s\n','Z levels');
for i = 1:zlevs
    fprintf(fid,'%i %f\n',[i -zdepth(i)]);
end
fprintf(fid,'%s\n','S levels');
fprintf(fid,'%f %f %f\n',[contrans pos constr]);
for i = 1:slevs
    fprintf(fid,'%i %f\n',[i sigma(i)]);
end
fclose(fid);

out = exist('vgrid.in','file');