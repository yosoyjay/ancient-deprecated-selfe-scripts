#!/usr/bin/python
#
# This script goes through a particle.bp built up from a backtrack run,
# removes those particles above the free surface, and adjusts the number of particles.
#
# Run folder names expected to adhere to my convention:  d-m-y-n_f
# Input: number of days in run (n)
# Output: particle.bp with particles above free surface removed
#
# jlopez 11/04/2010


import sys
import os
import shutil
#import subprocess
import commands

# Grab number of days in the run
if len(sys.argv) != 2:
	sys.exit("Usage: %s [number of days in the run n]" % sys.argv[0])
else:
	num_days = int(sys.argv[1])

# Run through particle.pth file
# 1. Throw out first line
# 2. Capture second line which gives the number of time steps
# 3. Capture third line with time and number of particles (should be 5000)
# 4. For i=1:num_time_steps
# 5.  For i=1:num_particles
# 6.    Read and save particle info to data_structure_n save to file particle_n.pth
# 7. end

pth_source_fid = os.open(path_to_pth, 'r');
null = pth_source_fid.readline()
num_time_steps = int(pth_source_fid.readline())

for i in range(1, num_time_steps+1):
	time_step = pth_source_fid.howeveryoureadonenumber
	parts = pth_source_fid.howeveryoureadonenumber

	# Preallocate space for array of particles 
	parts_array[3*parts]
	for j in range(1, parts):
		parts_array(j) = pth_source_fid.readline()

	file_name = 'particle_%i' % time_step
	output_fid = os.fopen(file_name, 'w')



for day in range(1,num_days+1):
	# 8 should be a varible for days in run + 1?
	run_length = num_days + 1 - day	
	path='7-%d-2010_%d_f' % (day, run_length)
	os.chdir(path)
	os.system('tail -n 5000 particle.bp | awk \'$5 >= 0\' > above_free_surface.dat')
	os.system('tail -n 5000 particle.bp | awk \'$5 < 0\' > new_particle.tmp')
	os.system('head -n 7 particle.bp > new_head.tmp')
	os.system('cat new_head.tmp new_particle.tmp > particle.bp.new')
	# If python is ever updated I can used the non-deprecated version
	#proc = subprocess.Popen('awk \'END{print NR}\' above_free_surface.dat', shell=True, stdout=subprocess.PIPE, )
	#num_afs = proc.communicate()[0]	
	num_afs = int(commands.getoutput('awk \'END{print NR}\' above_free_surface.dat'))
	print '%d particles above the free surface to be deleted from particle.bp in 7-%d-2010_%d_f' % (num_afs, day, run_length)
	new_part_num = 5000-num_afs
	os.system('sed -e "7s/5000/%d/" particle.bp.new > particle.bp' % new_part_num)
	print '%d particles now in particle.bp' % new_part_num
	os.system('rm *.tmp; rm *.new')
	os.chdir('..')

	
		
	
	
	
