% This function reads particles that ended up above the free surface
% at the end of a backtracking run ended up from parts_afs.dat and 
% uses those particle numbers to remove those particles from a copy
% of particle.bp. Returns a struct with a time of 0, x, y, z to adhere
% to my own imposed data structure standarization policy.
%
% jlopez 11/08/2010
%

function parts_minus_afs = remove_afs_bp(path_to_afs, path_to_bp)

afs_fd = fopen(path_to_afs);
if afs_fd<0
    disp('Error opening file listing particles above free surface');
    fg = [];
    return;
end

bp_fd = fopen(path_to_bp);
if bp_fd<0
    disp('Error opening particle.bp');
    fg = [];
    return;
end

% Read in all the info about the particles above free surface 
afs = load(path_to_afs);
afs_num = size(afs,1);
	
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


% Now delete the appropriate particles from parts
% Using A(n,1) = [] shifts all the values up, needing to account
% for that in the erase loop results in the -(i-1) term.
for i = 1:afs_num
    shift = i - 1;
	parts(afs(i,1)- shift,:) = [];
end

parts_minus_afs = struct( 'time', 0, 'x', parts(:,3), 'y', parts(:,4), ...
						  'z', parts(:,5) );


% Tidy up
fclose(afs_fd);
fclose(bp_fd);
