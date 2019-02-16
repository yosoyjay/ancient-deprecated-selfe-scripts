#!/usr/bin/python
""" This script is just a wrapper for readssh2b to generate Z0.out

Works for hi-res(15), 22, 26, and hi-res-26(27)
"""

import sys
import os
import datetime
import shutil
import subprocess

# Constants...
BIN_PATH = "/home/workspace/project/lopezj/bin/selfe_setup/"

def genZ0( startDate, grid ):

# Generate Z0.out for bctides.in.
# This version of readssh2b uses ssh.in for input instead of date.in.
# They are the same thing.
	dateIn = open("./ssh.in", "w")
	dateIn.write("%i %i %i 0\n" % (startDate[0], startDate[1], startDate[2]))    
	dateIn.write("0\n")
	if grid == 22 or grid == 26:
		dateIn.write("2 1 2\n")			# Number of boundaries that need Z0 
	else:
		dateIn.write("1 1\n")
	dateIn.close()

	try:
		_cmd_args = "%s/readssh2b" % BIN_PATH
		subprocess.check_call(_cmd_args)
	except CalledProcessError:
		print("Error: problem generationg Z0.out from readssh2b")
		raise
	
if __name__ == '__main__':
# Grab start and number of days 
	if len(sys.argv) != 3:
		sys.exit("Usage: %s [begin_date mm-dd-yyyy] [grid (15, 22, or 26)]" 
				 % sys.argv[0])
	else:
		start = sys.argv[1]
		grid  = int(sys.argv[2])
	
	if grid != 15 and grid != 22 and grid != 26:
		sys.exit("Grid must be 15, 22, or 26\n")
	 
	(month, day, year) = start.split('-')
	startDay   = int(day)
	startMonth = int(month)
	startYear  = int(year) 

# Python style dates for script.
	startDate = datetime.date(startYear, startMonth, startDay)
	startDate = startDate.timetuple()

	genZ0(startDate, grid)
