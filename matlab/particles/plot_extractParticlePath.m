function [ particles ] = extractParticlePath(  particlePath,   ...
												   timeDirection,  ...
												   timeSteps,      ...
												   particleNo )

% This function reads the particle.pth file and returns the relavent particle paths. 
% ! Only built to deal with a single day run.  timeSteps <= 96
%
% particlePath = the path to the particle.pth file
% timeDirection = [f if going forward] [b if going backward]
% timeSteps = the number of time steps in the backtracking run
% particleNo = the specific particle number
%
% jlopez 2/07/2011
% 

if(timeDirection ~= 'b' && timeDirection ~= 'f')
	error('plot_extractParticlePath:timeDirection', 'timeDirection must be b or f');
end

allParticlesFromPth = read_pth( particlePath );
%bp = read_bp ( ftBpPath );
% Need to offset by one to keep timesteps aligned with data

% extract the data of a particular particle from .pth and place in an array
index = 1;
if( timeDirection == 'b')
	% 97-timesteps*4 -> i.e. 97-1*4 = 93,94,95,96 = 4 positions -> 4 positions = 1 timesteps
	% there are 96 steps in the particle tracking.
	% need to preallocate if doing large amounts of particles.
	for i = (97-timeSteps):96
		particles(index,1) = allParticlesFromPth(i,1).x(particleNo,1);
		particles(index,2) = allParticlesFromPth(i,1).y(particleNo,1);
		particles(index,3) = allParticlesFromPth(i,1).z(particleNo,1);
		index = index + 1;
	end
else 
	for i = 1:timeSteps
		particles(index,1) = allParticlesFromPth(i,1).x(particleNo,1);
		particles(index,2) = allParticlesFromPth(i,1).y(particleNo,1);
		particles(index,3) = allParticlesFromPth(i,1).z(particleNo,1);
		index = index + 1;
	end
end

