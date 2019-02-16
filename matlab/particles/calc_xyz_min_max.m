function [ x_min x_max y_min y_max z_min z_max ] = calc_xyz_min_max( bt_part_path, ft_part_path, time_steps, part_num )

% This function take two sets of particles as arguments and returns an array
% of the min and max in the xyz directions.
%
% I used it for plotting.
%
% jlopez 12/28/2010
%

part_set_1 = read_pth( bt_part_path );
part_set_2 = read_pth( ft_part_path );

% Give reasonable default values
x_min = part_set_1(1,1).x(part_num,1); 
x_max = part_set_1(1,1).x(part_num,1);
y_min = part_set_1(1,1).y(part_num,1);
y_max = part_set_1(1,1).y(part_num,1);
z_min = part_set_1(1,1).z(part_num,1);
z_max = part_set_1(1,1).z(part_num,1);

% for centering the graph and setting an appropriate window.
% find x,y,z min and max
for i=1:time_steps
	x_min_temp = min([ part_set_1(i,1).x(part_num,1) part_set_2(i,1).x(part_num,1) ]);
	x_min = min([x_min_temp x_min]);
	
	x_max_temp = max([ part_set_1(i,1).x(part_num,1) part_set_2(i,1).x(part_num,1) ]);
	x_max = max([x_max_temp x_max]);

	y_min_temp = min([ part_set_1(i,1).y(part_num,1) part_set_2(i,1).y(part_num,1) ]);
	y_min = min([y_min_temp y_min]);
	
	y_max_temp = max([ part_set_1(i,1).y(part_num,1) part_set_2(i,1).y(part_num,1) ]);
	y_max = max([y_max_temp y_max]);
	
	z_min_temp = min([ part_set_1(i,1).z(part_num,1) part_set_2(i,1).z(part_num,1) ]);
	z_min = min([z_min_temp z_min]);
	
	z_max_temp = max([ part_set_1(i,1).z(part_num,1) part_set_2(i,1).z(part_num,1) ]);
	z_max = max([z_max_temp z_max]);
end	


x_min = x_min - x_min*.0005;	
x_max = x_max + x_max*.0005;
y_min = y_min - y_min*.0001;
y_max = y_max + y_max*.0001;
z_min = z_min - z_min*.0005;
z_max = z_max + z_max*.0005;


