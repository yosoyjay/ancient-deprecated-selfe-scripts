% This file opens a particle.bp file and stores the particle information
% in a struct with a time, x, y, z to adhere to my own imposed data 
% structure standardization policy.
%
% jlopez 11/16/2010
%

function particles_from_bp = read_bp(path_to_bp)

bp_fd = fopen(path_to_bp);
if bp_fd < 0
	disp('Error opening particle.bp.  Check the path');
	return;
end

% Read particle starting information from particle.bp and store
% in a array.
null_dev = fgetl(bp_fd);		% Throw away most header stuff
[null_dev, count] = fscanf(bp_fd,'%i',1);	
[null_dev, count] = fscanf(bp_fd,'%i',1);
[null_dev, count] = fscanf(bp_fd,'%i',1);
[null_dev, count] = fscanf(bp_fd,'%i %f %f', 3);
[null_dev, count] = fscanf(bp_fd,'%f %i %i %i %i %i', 6);
[parts_num, count] = fscanf(bp_fd, '%i', 1);		% The # of particles are useful	


% Read in the initial particle information
parts = zeros(parts_num,5);
for i = 1:parts_num
	temp = fscanf(bp_fd, '%i %i %f %f %f', 5);
	parts(i,1) = temp(1);
	parts(i,2) = temp(2);
	parts(i,3) = temp(3);
	parts(i,4) = temp(4);
    parts(i,5) = temp(5);
end

% Shove data into struct
particles_from_bp = struct( 'time', 0, 'x', parts(:,3), 'y', parts(:,4), ...
							'z', parts(:,5) );

% Clean up shop
fclose(bp_fd);

