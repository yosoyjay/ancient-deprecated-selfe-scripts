% Attept to draw regions for histogram data analysis

cd '/Users/jesse/Documents/selfe/practice/';
fg = hgrid2fg('hgrid.gr3');

% roughly define a box 
x_min=3.31e5;
x_max=3.68e5;
y_min=2.71e5;
y_max=3.09e5;

% define a region - find nodes within the region
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

pt_bndx=[E_bnd_x; NaN; isle_1_x; NaN; isle_2_x; NaN; isle_3_x; NaN; isle_4_x; NaN; isle_5_x; NaN; isle_6_x; NaN; isle_7_x; NaN; isle_8_x]; 
pt_bndy=[E_bnd_y; NaN; isle_1_y; NaN; isle_2_y; NaN; isle_3_y; NaN; isle_4_y; NaN; isle_5_y; NaN; isle_6_y; NaN; isle_7_y; NaN; isle_8_y];

%plot the region you've defined and compare to the boundary to make sure you aren't off by a point anywhere
clf;
plotbnd(fg);
outline = line(pt_bndx,pt_bndy);
set(gca,'xlim',[3.3e5 3.85e5],'ylim',[2.7e5 3.1e5]);


