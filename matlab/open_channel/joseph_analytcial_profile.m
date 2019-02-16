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
title('Analytic Vel. profile');
%ylim([0 H/10]);
subplot(2,1,2);
plot(niu,zp);
title('Analytic Viscosity profile');

%Output vvd.dat
fid=fopen('vvd.dat','w');
fprintf(fid,'%d\n',nvrt);
for i=1:nvrt
  fprintf(fid,'%d %f %f\n',i,niu(i),0.);
end
fclose(fid);

%Comp. with model
if(1==1)
  fid=fopen('../Scripts/station_scaled.bp','r');
  fgetl(fid);
  nsta=str2num(fgetl(fid));
  for i=1:nsta
    nums=str2num(fgetl(fid));
    xy(i,1:2)=nums(2:3);
  end
  fclose(fid);

%Load model
  m1_vel=load('../RUN14/vel_prof.dat'); %(1:nvrt2,dis(1:nvrt2),u(1:nvrt2),v(1:,nvrt2)) for all pts; dis is distance from bottom
  nvrt1=ceil(size(m1_vel,1)/nsta)
  m2_vel=load('../RUN13/vel_prof.dat');
  nvrt2=ceil(size(m2_vel,1)/nsta)
  m3_vel=load('../RUN15/vel_prof.dat');
  nvrt3=ceil(size(m3_vel,1)/nsta)

  figure(3);
  for i=2:nsta-1
    j1=(i-1)*nvrt1+1; %starting row for model results
    j2=(i-1)*nvrt2+1;
    j3=(i-1)*nvrt3+1;
    subplot(3,1,i-1); hold on;
    plot(u,zp,'r',m1_vel(j1:j1+nvrt1-1,3),m1_vel(j1:j1+nvrt1-1,2),'b',...
    m2_vel(j2:j2+nvrt2-1,3),m2_vel(j2:j2+nvrt2-1,2),'g',m3_vel(j3:j3+nvrt3-1,3),m3_vel(j3:j3+nvrt3-1,2),'k');
    title(strcat('(x,y)=(',num2str(xy(i,1:2)),')'));
    ylim([0 H]);
    %ylim([0 H/10]);
    set(gcf,'Color',[1 1 1]);
    if(i==nsta-1)
      legend('Analytical', ...
	strcat('Nz =',num2str(nvrt1)), ...
	strcat('Nz =',num2str(nvrt2)), ...
	strcat('Nz =',num2str(nvrt3)), ...
	'location', 'best');
    end
  end %for
end %if
