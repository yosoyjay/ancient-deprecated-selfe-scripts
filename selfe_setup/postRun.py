#!/usr/local/bin/python
""" This script creates hotstart files for all the times in the model
and cleans up the output directory and run directory.

Specifics:
- Generates last day of run hotstart, or optionally, hotstart files for all 
  the days of the run
- Deletes no files, all files but those required to continue a run, or bare
  bones files describing the run parameters
- Moves output files to the run directory out of the outputs directory to 
  enable ease of use with standard post-processing scripts

lopezj - 6/28/2011
"""

import os
import sys
import shutil
from optparse import OptionParser

# constant stuff and imported functions
# TODO: Fix the constants crap
HOME_PATH    = "/home/workspace/project/lopezj/"
SCRIPT_PATH = '/home/workspace/project/lopezj/scripts/ss_dev/'

sys.path.append(SCRIPT_PATH)
sys.path.append("%sbin/pexpect/" % HOME_PATH)
import pexpect
from continueRun import getOldRunData, getNumProcsUsed, genNewHotstart

# Hard coded source file for sflux and line numbers in param.in
# TODO: Make a param.in parser
SFLUX_SRC = 1
N_DAYS_LN = 151
IHFSKIP_LN = 262
DT_LN = 154
NTRACERS_LN = 24

def genHotstartFile(runDayNumber, ntracers):
	"""
	Generates a new hotstart file of day number n named hotstart.in_Day_n
	Assumes this is executed in a /run directory and moves into /run/outputs
	during the course of execution.

	"""

# Get old data
	(oldStartDate, oldNDays, old_ihfskip, time_step, ntracers) = getOldRunData()
	n_procs = getNumProcsUsed()

# Determine which hotstart output to combine 
	ihfskip = old_ihfskip*runDayNumber

# If the partial files do not exists, print error and then continue along.
	if os.path.exists("outputs/%i_0000_hotstart" % ihfskip):
		pass
	else:
		print("Unable to create hotstart file for day %i because the partial"
		      " file does not exist" % runDayNumber)
		return

# Generate new hotstart.in
	os.chdir('outputs')
	print "Generating hotstart for day %i" % runDayNumber

	if os.path.exists('combine_hotstart'):
		pass
	else:
		os.system('ln -s %s/bin/combine_hotstart ./' % HOME_PATH)

	child = pexpect.spawn('combine_hotstart', timeout= 240)
	child.expect('\n')
	child.sendline('%i %i' % (n_procs, ntracers))
	child.expect('\n')
#	child.expect('Input iteration # of the hotstart file (before _00*_hotstart) :\r\n')
	child.sendline('%i' % ihfskip)
	child.expect('\n')
	child.expect('\n')
	child.expect('Done')
	child.close()

# Clean up and move new hotstart to run folder
	shutil.move('hotstart.in', 'hotstart.in_Day_%i' % runDayNumber)
	os.chdir('..')


def postRun(hotstart, clean):
	""" Main part of the script.  Generates hotstart files, 
	moves the output files from outputs to run, and then deletes the output
	folder.
	"""
	(startDate, nDays, ihfskip, dt, ntracers) = getOldRunData()

# generate hostart file(s)
	if hotstart == 0:
		pass
	elif hotstart == 1:
		genHotstartFile(nDays, ntracers)
		os.system("mv outputs/hotstart.in* ./")
	elif hotstart == 2:
		for day in range(1,nDays+1):
			genHotstartFile(day, ntracers)
		os.system("mv outputs/hotstart.in* ./")
	else:
		sys.exit("usage: %s [hotstart option] [clean up option]\n"
				 "hotstart option:\n"
				 "\t0 - do not create hotstart.in\n"
				 "\t1 - create hotstart.in for last day of run\n"
				 "\t2 - create hotstart.in for every day of run\n"
				 "cleanup option:\n"
				 "\t0 - do not delete anything\n"
				 "\t1 - leave the files required to continue a run\n"
				 "\t2 - only leave bctides.in, param.in, *.ic, vgrid.in, and hgrid.in")

# Move output files to run directory. Lame. 
	os.system("mv outputs/?_[a-z]* ./")
	if nDays > 9:
		os.system("mv outputs/??_[a-z]* ./")

# Delete the goods
	if clean == 0:
		pass
	elif clean == 1:
		os.system("rm -f *_nu.in")
		os.system("rm -f fort.*")
		os.system("rm -f *.dat")
		os.system("rm -f outputs/*")
		os.system("rm -f *.o*")
		os.system("rm -f *.po*")
		os.system("rm -f autocombine.pl")
	elif clean == 2:
		os.system("rm -f *_nu.in")
		os.system("rm -f fort.*")
		os.system("rm -f *.dat")
		os.system("mv hgrid.gr3 hgrid~")
		os.system("rm -f *.gr3")
		os.system("mv hgrid~ hgrid.gr3")
		os.system("rm -f *.ll")
		os.system("rm -f *.o*")
		os.system("rm -f *.po*")
		os.system("rm -f autocombine.pl")
		os.system("rm -rf sflux")
		os.system("rm -rf outputs")
	else:
		sys.exit("usage: %s [hotstart option] [clean up option]\n"
				 "hotstart option:\n"
				 "\t0 - do not create hotstart.in\n"
				 "\t1 - create hotstart.in for last day of run\n"
				 "\t2 - create hotstart.in for every day of run\n"
				 "cleanup option:\n"
				 "\t0 - do not delete anything\n"
				 "\t1 - leave the files required to continue a run\n"
				 "\t2 - only leave bctides.in, param.in, *.ic, vgrid.in, *.th, and hgrid.gr3"
				 % sys.argv[0])
		
if __name__ == "__main__":
	if len(sys.argv) != 3:
		sys.exit("Usage: %s [hotstart option] [clean up option]\n"
				 "Hotstart option:\n"
				 "\t0 - Do not create hotstart.in\n"
				 "\t1 - Create hotstart.in for last day of run\n"
				 "\t2 - Create hotstart.in for every day of run\n"
				 "Cleanup option:\n"
				 "\t0 - Do not delete anything\n"
				 "\t1 - Leave the files required to continue a run\n"
				 "\t2 - Only leave bctides.in, param.in, *.ic, vgrid.in, and hgrid.in"
				 % sys.argv[0])
	else:
		hotstart = int(sys.argv[1])
		clean    = int(sys.argv[2])

	if(hotstart<0 or hotstart>2 or clean<0 or clean>2):
		sys.exit("Usage: %s [hotstart option] [clean up option]\n"
				 "Hotstart option:\n"
				 "\t0 - Do not create hotstart.in\n"
				 "\t1 - Create hotstart.in for last day of run\n"
				 "\t2 - Create hotstart.in for every day of run\n"
				 "Cleanup option:\n"
				 "\t0 - Do not delete anything\n"
				 "\t1 - Leave the files required to continue a run\n"
				 "\t2 - Only leave bctides.in, param.in, *.ic, vgrid.in, and hgrid.in" %
				 sys.arv[0])

	postRun(hotstart, clean)
