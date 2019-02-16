function vgrid=read_vgrid(file);

data=char(textread([file],'%s','delimiter','\n'));
s=str2num(data(1,:));
vgrid.nvrt=s(1);
vgrid.kz=s(2);
vgrid.h_s=s(3);

I=1;
for k=3:2+vgrid.kz
    s=str2num(data(k,:));
    vgrid.Zlevelindex(I)=s(1);
    vgrid.Zcoord(I)=s(2);
    I=I+1;
end

k=2+vgrid.kz+1;
s=str2num(data(k+1,:));

vgrid.hc=s(1);
vgrid.theta_b=s(2);
vgrid.theta_f=s(3);

k=k+2;
I=1;
for r=k:k+vgrid.nvrt-vgrid.kz
    s=str2num(data(r,:));
    vgrid.Slevelindex(I)=s(1);
    vgrid.Scoord(I)=s(2);
     I=I+1;
end