%Add m-elio to the matlab path

datafile = 'habitat_sal_AVG.61';
addpath '/usr/local/cmop/matlab/cmop/m-elio';
gr.hgrid=gr_readHGrid('hgrid.gr3');
h=sz_readHeader(datafile);
it = 1;
[d ts]=sz_readTimeStep(h,it);
axis equal;
xlim([309491.906309 406390.731852]);
ylim([249063.038661 315713.024484]);
gr_plot(gr.hgrid,d);
%quit;
