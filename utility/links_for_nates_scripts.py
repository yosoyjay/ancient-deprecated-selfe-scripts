#!/usr/local/bin/python
""" This scripts creates a directory of links in the form that Nate's
scripts expect.  

/top_level  <-- Run here - A new directory is created here
	/run 	
		/outputs
		/hgrid.gr3
		...
	/grid   <-- Links are created here. grid == 22, 26, etc... 
		/yyyy-ddd/
			/1_salt.63
			/?_blah.6?
			...

After this is run, you can use Nate's scripts.  For salinity cross section,
I am running scripts from top_level/ with that as base directory.  I copied
a shape directory from Nate's home directory.

NOTE! You must enter the correct starting date of the run because this is 
how the dates are calculated.  

Grid number is not important, it just ends up as the top of the top level
of the symlinks created for Nate's stuff.

jlopez - 05/09/2011
"""
#TODO
#1. Fix call to create symlinks to use call check
#2. Cleanup names
#3. Pull the date out of bctides.in
#4. Pull the number of days out of param.in


import sys
import datetime
import os
import glob

def create_nates_links(startDate, nDays, grid):
	top_dir = grid
	os.mkdir(top_dir)
	for i in range(nDays):
		run_day = "%s-%s" % (startDate[0], startDate[7]+i) 
		os.mkdir(top_dir + "/" + run_day)
		os.mkdir(top_dir + "/" + run_day + "/run")
		os.chdir(top_dir + "/" + run_day + "/run")
		os.symlink("../../../run/hgrid.gr3", "./hgrid.gr3")
		if os.path.exists("outputs"):
			os.system("ln -s ../../../run/outputs/%i_[a-z]* ./" % (i+1))
		else:
			os.system("ln -s ../../../run/%i_[a-z]* ./" % (i+1))

# Rename files if they are of type n_* to 1_* as that is what the scripts expect
		if i > 0:
			for file_name in glob.iglob("%i_*" % (i+1)):
				day_number, rest_of_name = file_name.split('_')
				os.rename(file_name, "1_%s" % rest_of_name) 
		os.chdir("../../../")

if __name__ == "__main__":
	if len(sys.argv) != 4:
		sys.exit("Usage: %s [start date mm-dd-yyyy] [number of days n]"
		         " [grid number (22,26,etc.)]" % sys.argv[0])

	start = sys.argv[1]
	nDays = int(sys.argv[2])
	grid  = sys.argv[3]

# deal with start and ending days
	(month, day, year) = start.split('-')
	startDay   = int(day)
	startMonth = int(month)
	startYear  = int(year) 

# Python style dates for script.
	startDate = datetime.date(startYear, startMonth, startDay)
	startDate = startDate.timetuple()

	create_nates_links(startDate, nDays, grid)
