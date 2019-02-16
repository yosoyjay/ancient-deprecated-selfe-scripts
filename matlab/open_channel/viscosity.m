%Calculated analytical values, and compare with model
%Input vgrid: no comments please; assume pure sigma 
clear all; close all;

H=20; %total depth
%Read in vgrid
fid=fopen('vgrid.in','r');
nums=str2num(fgetl(fid));
if(nums(2) ~=1); error('No SZ pls'); end;
nvrt=nums(1);
fgetl(fid); fgetl(fid); fgetl(fid);
nums=str2num(fgetl(fid));
hc=nums(1); theta_b=nums(2); theta_f=nums(3);
if(theta_f > 1.e-4); error('Pure sigma pls'); end;
for i=1:nvrt
  nums=str2num(fgetl(fid));
  sig(i)=nums(2);
  zp(i)=H*(1+sig(i));
end
fclose(fid);
figure(1);
plot(ones(nvrt,1),zp,'ko');

%Const.
gam=2e-4; %slope
Cd=1.e-2; %for 2D (Manning)
z1=1.e-3; %small const.
grav=9.81;
dahv=sqrt(grav*gam*H/Cd); %average vel
lam=grav*gam/dahv*(log(H/z1)-1);
format long;
disp('dahv lam='); [dahv lam]

%zp=0:H/1000:H; %distance from bottom
u=grav*gam/lam*log(zp/z1+1);
niu=lam*(H-zp).*(zp+z1);
Cd3D=grav*gam*H/u(2)/u(2); %3D drag
%Compute aver. vel. from u
dahv2=trapz(zp,u)/H;
disp('min/max visco., visco. at bottom='); [min(niu) max(niu) niu(1)]
disp('3D Cd='); Cd3D
disp('Depth-ave. vel from 3D profiel='); dahv2

figure(2);
subplot(2,1,1);
plot(u,zp,'k-o');
title('Vel. profile');
%ylim([0 H/10]);
subplot(2,1,2);
plot(niu,zp);
title('Viscosity profile');
