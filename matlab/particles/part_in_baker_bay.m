% This function takes particle positions and returns the indices of the
% of the particles in each respective region.
%
% +bb = Baker Bay     
% +yb = Youngs Bay
% ?cb = Calathmet Bay
% ?gb = Greys Bay
% ?em = Estuary Mouth
% nc = North Channel
% sc = South Channel
% ws = Washington Shelf
% os = Oregon Shelf
% do = Deep Ocean (Off the continental shelf)
% 
% TODO:
% 1. Finish creating boundaries for each region

function [bb_parts, yb_parts, cb_parts, gb_parts,...
		  em_parts, nc_parts, sc_parts, ws_parts,...
		  os_parts, do_parts] = part_in_baker_bay(particles);

%
% Region: Baker Bay
%

% Boundaries
x_1_min = 3.359e5
x_1_max = 3.438e5
pt_nw = [3.37e5 2.948e5]
pt_se = [3.459e5 2.924e5]

% Load particles into x and y vectors here
x_parts = particles.x;
y_parts = particles.y;

% Filter particles to appropriate variable
bb_parts = find(x_parts > x_1_min & ... 
				x_parts < x_1_max & ...
				y_parts > linear_func(x_parts, pt_nw, pt_se))

%
% Region: Young's Bay
%

% Boundaries
x_min = 3.453e5
x_max = 3.542e5
y_max = 2.860e5

% Filter particles
yb_parts = find(x_parts > x_min & ... 
				x_parts < x_max & ...
				y_parts < y_max)

%
% Region: Calathmet Bay
%

% Boundaries
x_min = 3.578e5;
y_max = 2.881e5;

% Filter particles
cb_parts = find(x_parts > x_min & ... 
				y_parts < y_max);

%
% Region: Grays Bay
%

% Boundaries
x_min = 3.586e5;
x_max = 3.648e5;
y_min = 2.881e5;

% Filter particles
gb_parts = find(x_parts > x_min & ... 
				x_parts < x_max & ...
				y_parts > y_min)

%
% Region: Estuary Mouth 
%

% Boundaries
x_min = 3.342e5; 
x_max = 3.437e5;

% Filter particles
em_parts = find(x_parts > x_min & ... 
				x_parts < x_max )


%
% Region: North Channel
%

% Boundaries
x_1_min = 3.437e5;
x_1_max = 0;
y_1_min = 0;
y_1_max = 0;
x_2_min = 0;
x_2_max = 0;
y_2_min = 0;
y_2_max = 0;
x_3_min = 0;
x_3_max = 0;
y_3_min = 0;
y_3_max = 0;

% Filter particles
%nc_parts = find(x_parts > x_min & ... 
%				x_parts < x_max )
nc_parts = 0;

%
% Region: South Channel
%

% Boundaries
x_1_min = 3.437e5;
x_1_max = 0;
y_1_min = 0;
y_1_max = 0;
x_2_min = 0;
x_2_max = 0;
y_2_min = 0;
y_2_max = 0;
x_3_min = 0;
x_3_max = 0;
y_3_min = 0;
y_3_max = 0;

% Filter particles
%sc_parts = find(x_parts > x_min & ... 
%				x_parts < x_max )
sc_parts = 0;
%
% Change to ocean part of the grid.
%

%
% Region: Washington Shelf
%

% Boundaries
x_min = 3.437e5;
%x_max = DOMAIN_BOUNDARY;

% Filter particles
%ws_parts = find(x_parts > x_min & ... 
%				x_parts < x_max )
ws_parts = 0;

%
% Region: Oregon Shelf
%

% Boundaries
x_1_min = 0;
x_1_max = 3.438e5;
y_1_min = 2.6e5
y_1_max = 2.806e5;
x_2_min = 0;
x_2_max = 3.421e5;
y_2_min = 2.806e5;
y_2_max = 2.883e5;
x_3_min = 0;
x_3_max = 3.396e5;
y_3_min = 2.883e5;
y_3_max = 2.909e5;

% Filter particles
%os_parts = find(x_parts > x_min & ... 
%				x_parts < x_max )
os_parts = 0;


%
% Region: Deep Ocean (Off the shelf) 
%

% Boundaries
x_min = 3.437e5;
%x_max = DOMAIN_BOUNDARY;

% Filter particles
%do_parts = find(x_parts > x_min & ... 
%				x_parts < x_max )
do_parts = 0;
