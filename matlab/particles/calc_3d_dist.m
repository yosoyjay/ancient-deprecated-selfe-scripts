% Inputs: - Two points in R3
% Ouput: Euclidean distance between two points

function [distance] = calc_3d_dist( x1, x2, y1, y2, z1, z2 )

diff_x = abs(x1 - x2);
diff_y = abs(y1 - y2);
diff_z = abs(z1 - z2);

distance = sqrt( diff_x*diff_x + diff_y*diff_y + diff_z*diff_z );
