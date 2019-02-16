% This function compares the particle positions after going through
% backtracking and then forward tracking. Using the distance formula
% it calculates the distance between each particles origin and 
% destiniation position returned in a vector.
% 
% The particle.bp file needs to be from the beginning of the backtracking run.
% The particle.pth file needs to be from the forward tracking results.
%
% Depends on:
% read_pth.m   - Reads in particle.pth files
% remove_afs.m - Reads in original particle.bp from backtracking and removes the
%                particles that ended up above the free surface.
% calc_dist.m  - Calculates distance between each particle in horizontal plane
% 				 between two sets of particles.
%
% jlopez 11/09/2010
%

function diff_orig_dest = delta_bt_ft(path_to_particle_bp, path_to_particle_pth, ...
									  path_to_afs_dat)

% Read in the particles from their destination runs - read from particle.pth from
% forward tracking run.
temp_parts_destination = read_pth(path_to_particle_pth);

% Read in particles original position and remove the particles that ended up above 
% the free surface at the end of the backtracking run.
parts_origin = remove_afs(path_to_afs_dat, path_to_particle_bp);

% Calculate the difference between each particles origin and destination.
% Different extraction methods because the data structures are different for each
% source.  This just gets the last structure which holds the time and
% position of the particles at the end of the run.  
dest_index = size(temp_parts_destination,1)
parts_destination = temp_parts_destination(dest_index)

% Calculates the distance in the horizontal plane between two sets of particles
diff_orig_dest = calc_dist(parts_origin, parts_destination); 
