% This is a script to specify random positions for the particles in the region 
% under the influece of the ocean as determined by the extent of salinity in the 
% estuary from 1997-2010.  The positions are distributed with a uniform random
% distribution in x and Y  and then given a random depth (not uniformly
% distributed 3-dimensionally).
%
% Modified version of a script bequethed to me from Scott Durski

working_dir = '/Users/jesse/Documents/selfe/practice';

cd '/Users/jesse/Documents/selfe/practice';

fg = hgrid2fg('hgrid.gr3');

Num_pts=5000;
num_rel_days=1;
%ret_rate=0.1;   					% release particles at 1/8th of the locations per day
ret_rate=1.0;
rel_per_day=1;   					% release new particles rel_per_day times a day
num_new_drft=fix(ret_rate/rel_per_day*Num_pts);		% ? Fix rounds numbers towards 0 to make integers
tot_drft=(num_rel_days)*num_new_drft;		

% roughly define a box 
x_min=3.31e5;
x_max=3.68e5;
y_min=2.71e5;
y_max=3.09e5;

% define a region - find nodes within the region?
xDbnd=fg.x(fg.bnd(:,1));
yDbnd=fg.y(fg.bnd(:,1));
indx=find(xDbnd>x_min & xDbnd<x_max & yDbnd>y_min & yDbnd<y_max);

% isle 1
isle_1_x=[xDbnd(indx(649:689)); xDbnd(indx(649))];
isle_1_y=[yDbnd(indx(649:689)); yDbnd(indx(649))];

% isle 2
isle_2_x=[xDbnd(indx(690:714)); xDbnd(indx(690))];
isle_2_y=[yDbnd(indx(690:714)); yDbnd(indx(690))];

% isle 6
isle_6_x=[xDbnd(indx(606:616)); xDbnd(indx(606))]; 
isle_6_y=[yDbnd(indx(606:616)); yDbnd(indx(606))]; 

% isle 4
isle_4_x=[xDbnd(indx(617:648)); xDbnd(indx(617))];
isle_4_y=[yDbnd(indx(617:648)); yDbnd(indx(617))];

% isle 5
isle_5_x=[xDbnd(indx(715:758)); xDbnd(indx(715))];
isle_5_y=[yDbnd(indx(715:758)); yDbnd(indx(715))];

% isle 3
isle_3_x=[xDbnd(indx(759:776)); xDbnd(indx(759))];
isle_3_y=[yDbnd(indx(759:776)); yDbnd(indx(759))];

% isle 8
isle_8_x=[xDbnd(indx(777:796)); xDbnd(indx(777))];
isle_8_y=[yDbnd(indx(777:796)); yDbnd(indx(777))];

% isle 7
isle_7_x=[xDbnd(indx(568:605)); xDbnd(indx(568))];
isle_7_y=[yDbnd(indx(568:605)); yDbnd(indx(568))];

% estuary boundary
E_bnd_x=[xDbnd(indx(47:533)); xDbnd(indx(47))]; 
E_bnd_y=[yDbnd(indx(47:533)); yDbnd(indx(47))]; 

BB_bndx=[E_bnd_x; NaN; isle_1_x; NaN; isle_2_x; NaN; isle_3_x; NaN; isle_4_x; NaN; isle_5_x; NaN; isle_6_x; NaN; isle_7_x; NaN; isle_8_x]; 
BB_bndy=[E_bnd_y; NaN; isle_1_y; NaN; isle_2_y; NaN; isle_3_y; NaN; isle_4_y; NaN; isle_5_y; NaN; isle_6_y; NaN; isle_7_y; NaN; isle_8_y];

%plot the region you've defined and compare to the boundary to make sure you aren't off by a point anywhere
clf;plotbnd(fg);
hl2=line(BB_bndx,BB_bndy);
set(gca,'xlim',[3.3e5 3.85e5],'ylim',[2.7e5 3.1e5])

% generate a box that nicely encompasses this region...
xmax=max(BB_bndx);xmin=min(BB_bndx);xrng=xmax-xmin;
ymax=max(BB_bndy);ymin=min(BB_bndy);yrng=ymax-ymin;

% find all nodes within this box for later search
i_inBx=find(fg.x>xmin & fg.x<xmax & fg.y>ymin & fg.y < ymax);
% find all elements associated with these nodes for later triangulations
for i_nd=1:length(i_inBx),
    inBx(i_nd).elems=find(fg.e(:,1)==i_inBx(i_nd) | fg.e(:,2)==i_inBx(i_nd) | fg.e(:,3)==i_inBx(i_nd) ) ;
end


%reltime=ones([Num_pts 1])*0.0;
%zdep=zeros([Num_pts 1]);

% for the hires db14 grid we don't want too many drifters so take half of the
% ones located in bakers bay

%load /Users/sdurski/Selfe/Particle_tracking/SS_2009/MAT_data/Elev_peaks_Saturn_1

%for irs=1:length(ht_peak)
for irs=1:1,
%for irs=hires_set,
    
  %offset_rel=timet(ht_time(irs));  % time offset for first release of the day
  offset_rel = 0;
  
  %if irs<10 
  %  filename=sprintf('/Users/sdurski/Selfe/Particle_tracking/SS_2009/Ensembles/Passive/BP_As_files/particle.bp.r00%i_As',irs)
  %elseif irs<100
  %    filename=sprintf('/Users/sdurski/Selfe/Particle_tracking/SS_2009/Ensembles/Passive/BP_As_files/particle.bp.r0%i_As',irs)
  %else
  %    filename=sprintf('/Users/sdurski/Selfe/Particle_tracking/SS_2009/Ensembles/Passive/BP_As_files/particle.bp.r%i_As',irs)
  %end
  %fid=fopen(filename,'wt');
  fid=fopen('/Users/jesse/Documents/selfe/practice/particle.bp', 'wt');
 
 % Write BP file
  fprintf(fid,'Mouth of the river to the extent of maximum salinity\n');
  fprintf(fid,'%i\n',1);
  fprintf(fid,'%i\n',1);
  fprintf(fid,'%i\n',0);  	% stiffness
  fprintf(fid,'%i %f %f\n',1, -124, 46.25);
 %fprintf(fid,'%f %i %i %i %i %i\n',0.1,fix(timet(ht_time(irs)))+31,150,72,576,100);
  fprintf(fid,'%f %i %i %i %i %i\n',0.1,1,90,10,960,100);
  fprintf(fid,'%i\n',tot_drft);

  itot=0;
  for it=1:num_rel_days,
    for ipd=1:rel_per_day,
      reltime=((it-1)+(ipd-1)/rel_per_day+offset_rel)*86400;  
 % generate two long random vectors and extract the x,y
 % locations that are within the region.
      xrnds=rand([Num_pts*3 1])*xrng+xmin;
      yrnds=rand([Num_pts*3 1])*yrng+ymin;
      ipt_inreg=find(inpolygon(xrnds,yrnds,BB_bndx,BB_bndy)==1);
 % find the element each point is within,
 % find the water depth at that point in that element -
 % eliminate particles in dry elements
 % randomly choose Num_pts particle locations from the set.
      for ip=1:length(ipt_inreg),
 % from a subset of all nodes, find the closest to the particle location
         xpt(ip)=xrnds(ipt_inreg(ip));
         ypt(ip)=yrnds(ipt_inreg(ip));
         dist=sqrt((xpt(ip)-fg.x(i_inBx)).^2+(ypt(ip)-fg.y(i_inBx)).^2);
         ind_close=find(dist==min(dist));
 % search the elements that node is associated with and find the one that 
 % triangulates the particle location.   
 % I think the way to do this is to sum the areas of the three triangles
 % formed by the point and each node in each candidate element and find
 % the case where the sum equals the area of the element.
 % Heron's formula calculates the area of a triangle given the lengths of
 % the three sides.
 % loop over candidate elements to find which one the point is in...
        el_cand=inBx(ind_close).elems;
        diff_ar=[];
        for ie=1:length(el_cand)
          iel_cnd=el_cand(ie);
 % calculate the area of the three triangles formed between the nodes of 
 % this element and the point of interest
          nd(1)=fg.e(iel_cnd,1);nd(2)=fg.e(iel_cnd,2);nd(3)=fg.e(iel_cnd,3);
          for is=1:3
              sp(is)=sqrt((xpt(ip)-fg.x(nd(is))).^2+(ypt(ip)-fg.y(nd(is))).^2);
          end
          snd(1)=sqrt((fg.x(nd(1))-fg.x(nd(2))).^2+(fg.y(nd(1))-fg.y(nd(2))).^2);
          snd(2)=sqrt((fg.x(nd(2))-fg.x(nd(3))).^2+(fg.y(nd(2))-fg.y(nd(3))).^2);
          snd(3)=sqrt((fg.x(nd(3))-fg.x(nd(1))).^2+(fg.y(nd(3))-fg.y(nd(1))).^2);
 % calculate three semi-perimeters and total semiperimeter of the element 
          per(1)=(sp(1)+sp(2)+snd(1))*0.5;
          per(2)=(sp(2)+sp(3)+snd(2))*0.5;
          per(3)=(sp(3)+sp(1)+snd(3))*0.5;
          per_el=(snd(1)+snd(2)+snd(3))*0.5;
 % calculate areas using Heron's formula
          ar(1,ie)=sqrt(per(1)*(per(1)-sp(1)).*(per(1)-sp(2)).*(per(1)-snd(1)));
          ar(2,ie)=sqrt(per(2)*(per(2)-sp(2)).*(per(2)-sp(3)).*(per(2)-snd(2)));
          ar(3,ie)=sqrt(per(3)*(per(3)-sp(3)).*(per(3)-sp(1)).*(per(3)-snd(3)));
          ar_el(ie)=sqrt(per_el*(per_el-snd(1)).*(per_el-snd(2)).*(per_el-snd(3)));
          diff_ar(ie)=ar_el(ie)-(ar(1,ie)+ar(2,ie)+ar(3,ie));
        end
 % the minimum abs(difference) gives the element we want.
        imn=find(abs(diff_ar)==min(abs(diff_ar)));
        pt_elem=el_cand(imn);
        nd(1)=fg.e(pt_elem,1);nd(2)=fg.e(pt_elem,2);nd(3)=fg.e(pt_elem,3);
 % the water depth at the particle position can then be calculated as
        zdep=(fg.z(nd(1)).*ar(2,imn)+fg.z(nd(2)).*ar(3,imn)+  ... 
             fg.z(nd(3)).*ar(1,imn))./ar_el(imn);
        zfault(ip)=min(0,zdep);
 % and the random vertical position of the particle can be calculated.
        zpt(ip)=rand(1).*(zdep);
      end
      in_zf=find(zfault~=0);
 % for multiple releases, check about writing 'ip' below...
      for ip=1:Num_pts,
       fprintf(fid,'%i   %i   %f   %f   %f\n',ip,    ...
                   reltime,xpt(in_zf(ip)),ypt(in_zf(ip)),zpt(in_zf(ip)));
      end
    end
  end
  fclose(fid);
  irs;
end


