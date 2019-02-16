% This function reads in *.pth files and seperates the position for each 
% particle (x,y,z) by time.  This input to the function is the path to a
% particle file and the output is a structure of the particles positions.
% This function assumes you know a good estimate of how many particles 
% are in the run.
%
% jlopez 11/8/2010

function [particle_positions] = particle_pos(path_to_pth_file)

[fid, error] = fopen(path_to_pth_file);

if fid<0
	disp(error)
end

%Drogues blah is just trashed
[drogues_blah, count] = fscanf(fid,'%s',1);  
[time_steps, count] = fscanf(fid,'%i',1);
step = 0;

% Preallocate array space where particle information will reside to avoid
% many reallocations while I shove structs in the array.  There should be a
% time_step x 1 array returned.
z_val(1:5000) = 0;
fake_struct = struct( 'time', 10000, 'x', z_val, 'y', z_val, 'z', z_val );
particle_positions = repmat(fake_struct, time_steps, 1);

while count>0
	step = step+1;

	% Read in information about the time step [ time and number of particles ] 
	% and save it in time and part_num
	[step_info, count] = fscanf(fid,'%f %f',2);
	if count==0
		path.x = x;
		path.y = y;
		path.z = z;
		path.time = time;
		return;
	end
	time(step) = step_info(1);
	part_num = step_info(2);

	% Read in each particles information (part_num, x, y, z) and save it to a 
	% structure put it in the array and go away.
	for i = 1:part_num
		[part_info, count] = fscanf(fid,'%i %f %f %f', 4);	
		if count~=4
			break
		end
		x(i,1) = part_info(2);
		y(i,1) = part_info(3);
		z(i,1) = part_info(4);
	end

	s = struct( 'time', step_info(1), 'x', x, 'y', y, 'z', z); 
	particle_positions(step,1) = s;	

end

% Be sure to tidy up
fclose(fid);
