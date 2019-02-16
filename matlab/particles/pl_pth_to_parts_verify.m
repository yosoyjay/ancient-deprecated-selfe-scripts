function [ btParticles ftParticle ] = extractParticlePaths( btPthPath,   ...
           									 				ftPthPath,   ...
									  						ftBpPath,    ...
												  			timeSteps,   ...
												  			particleNo )

% This function reads the backtracking.pth file and the forwardtracking.pth file
% from a bt and ft run and returns the appropriate positions to make a comparision
% between the two runs.
%
% timeSteps = the number of time steps in the backtracking run
% particleNo = the specific particle number
%
% jlopez 12/29/2010
%

% Read in all the relevant data
allBtParticlesFromPth = read_pth( btPthPath );
allBtParticlesFromPth = read_pth( ftPthPath );
bp = read_bp ( ftBpPath );

index = 1;
% 97-timeSteps*4 -> i.e. 97-1*4 = 93,94,95,96 = 4 positions -> 4 positions = 1 timeSteps
% there are 96 steps in the particle tracking.
% need to preallocate if doing large amounts of particles.
for i = (97-timeSteps):96;
	btParticles(index,1) = allBtParticlesFromPth(i,1).x(particleNo,1);
	btParticles(index,2) = allBtParticlesFromPth(i,1).y(particleNo,1);
	btParticles(index,3) = allBtParticlesFromPth(i,1).z(particleNo,1);
	index = index + 1;
end

% Need to get starting position for particle.bp to compare to end position of
% the forward tracking.  This is quite silly, ptrack should print starting position.
% Then I ignore the last position because it overshoots the backtracking by one
% time step.
% Note:  The forward tracking run must start at 900 seconds to properly line up
% with bt.
ftParticle(1,1) = bp.x(particleNo,1);
ftParticle(1,2) = bp.y(particleNo,1);
ftParticle(1,3) = bp.z(particleNo,1);

index = 2;
for i = 2:timeSteps;
	ftParticle(index,1) = allBtParticlesFromPth(i,1).x(particleNo,1);
	ftParticle(index,2) = allBtParticlesFromPth(i,1).y(particleNo,1);
	ftParticle(index,3) = allBtParticlesFromPth(i,1).z(particleNo,1);
    index = index + 1;
end

