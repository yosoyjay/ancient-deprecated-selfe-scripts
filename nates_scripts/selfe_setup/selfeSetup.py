#!/usr/local/bin/python
"""This script sets up the folders and files for a SELFE run. 

The run is set up assuming hotstart and no tracers.  Its creates a baseline
of files and settings that are commonly used.  Any changes to be made must
completed manually after the set up.

Note:  Some modules are imported immediately before the functions are called
due to dependencies fulfilled only at that point and poor planning on my
part.
"""
import os
import sys
import shutil
import datetime 
import subprocess
from optparse import OptionParser

# constant stuff and imported functions
SCRIPT_PATH = '/home/workspace/project/lopezj/scripts/selfe_setup/'
sys.path.append(SCRIPT_PATH)
SFLUX_SRC = 1

def selfeSetup(startDate, n_days, grid, n_tracers, email=''):

	endDate   = startDate + datetime.timedelta(days=(n_days-1))
	startDate = startDate.timetuple()		
	endDate   = endDate.timetuple()

# Make directories
	print 'Creating run directories'
	try:
		run_dir = "%i-%s-%s_db%i_%s_days" % (startDate[0],                   
			str(startDate[1]).zfill(2), str(startDate[2]).zfill(2), grid,   
			str(n_days).zfill(2))
		os.mkdir("%s" % run_dir) 
		os.mkdir("%s/run" % run_dir) 
		os.mkdir("%s/run/outputs" % run_dir) 
		os.mkdir("%s/post" % run_dir) 
		os.chdir("%s/run" % run_dir) 
	except OSError:
		print("Error: problem creating run directories")
		raise

# copy *.gr3 files from canonical versions
# Grid 27 is changed to grid = 26 here because it has the same boundary conditions
# and that is all we are concerned about for the rest of the script.
	print '\nlinking to gr3 files for db%i' % grid
	if(grid == 15):
		ret_code = os.system('ln -s %s/hi_res/* ./' % SCRIPT_PATH)
		if ret_code != 0:
			print("Error: problem creating symlinks of grid files to run directory")
			raise IOError
	elif(grid == 22):
		#os.symlink('%s/db22/*' % SCRIPT_PATH, './')
		ret_code = os.system('ln -s %s/db22/* ./' % SCRIPT_PATH)
		if ret_code != 0:
			print("Error: problem creating symlinks of grid files to run directory")
			raise IOError
	elif(grid == 26):
		#os.symlink('%s/db26/*' % SCRIPT_PATH, './')
		ret_code = os.system('ln -s %s/db26/* ./' % SCRIPT_PATH)
		if ret_code != 0:
			print("Error: problem creating symlinks of grid files to run directory")
			raise IOError
	else: 
# High resolution version of 26, set up grid stuff then the rest is based on 26.
		#os.symlink('%s/db22/*' % SCRIPT_PATH, './')
		ret_code = os.system('ln -s %s/db26_div4/* ./' % SCRIPT_PATH)
		if ret_code != 0:
			print("Error: problem creating symlinks of grid files to run directory")
			raise IOError
		grid = 26

# generate sflux - using nam (1) as a source - whatever that means...
	print '\ngenerating sflux\n'
#	try:
#		_cmd_args = ["%s/make_sflux_links.csh" % SCRIPT_PATH,
#			str(startDate[0]), str(startDate[1]), str(startDate[2]), 
#			str(endDate[0]), str(endDate[1]), str(endDate[2])]
#		subprocess.check_call(_cmd_args)
#	except:
#		print("Error: Problem generating sflux.")
#		raise
	command = "%s/make_sflux_links.csh" % SCRIPT_PATH
	return_code = subprocess.call([command, str(SFLUX_SRC), 
		str(startDate[0]), str(startDate[1]), str(startDate[2]),
		str(endDate[0]), str(endDate[1]), str(endDate[2])])
	if return_code != 0:
		sys.exit("Error: Problem generating sflux.")
			   
# generate hotstart.in assuming ocean is only boundary that is of flux type 4.
# first create date.in file
	try:
		dateIn = open("date.in", "w")
		try:
			if(grid==15):
				dateIn.write('1 1 1\n')		# 1 open boundary of type 4
			else:
				dateIn.write('1 1 2\n')		# 2 open boundary of type 4?
			dateIn.write('%i %i %i 0\n' % 
				(startDate[0], startDate[1], startDate[2]))    
			dateIn.write('0\n')			 	# ireduce
			dateIn.write('SELFE\n')			# selfe or elcirc
			dateIn.write('1.e6 33.\n')		# ht(depth), salinitymin 
			dateIn.write('%s\n' % n_tracers)# number of tracers
		finally:
			dateIn.close()
	except IOError:
		print("Error: unable to open date.in for writing\n")
		raise

# call to get temp.th and flux.th because temp.th is needed for hotstart.in
	print '\nGenertating flux.th and temp.th\n'
	from genFluxTemp  import genFluxTemp 
	genFluxTemp(startDate, endDate, grid)

# call pragram to generate hotstart.in from ncom
# Need to use readncom8b and not a because a doesn't work.
	print '\nGenerating hotstart.in\n'
	try:
		_cmd_args = "%s/bin/readncom8b" % SCRIPT_PATH
		subprocess.check_call(_cmd_args)
	except:
		print("Error: problem generating hotstart.in")
		raise
#	retrun_code = subprocess.call('%s/bin/readncom8b' % SCRIPT_PATH)
#	if return_code != 0:
#		sys.exit("Error: problem generating hotstart.in")

# Change date.in and call readncom to generate nudging files from NCOM
	print '\nGenerating nudging files\n'
	try:
		_cmd_args = ["sed", "-i", "2s/.$/%i/" % n_days, "date.in"]  
		subprocess.check_call(_cmd_args)
	except:
		print("Error: Problem editing date.in for nudging files")
		raise

	print '\nGenerating tidal data\n'
	from genTidalData import genTidalData
	genTidalData(startDate, grid)

	print 'Generating Z0 values\n'
	from genZ0 import genZ0
	genZ0(startDate, grid)

	print '\nGenerating bctides.in\n'
	from genBCTides   import genBCTides
	genBCTides(startDate, grid)

# Copy param.in and change number of days for run and turn on writing tracer 
# output files.  Only turns on two.
# Copy pelfe binary and job submit script.
	try:
		shutil.copy('%s/param.in' % SCRIPT_PATH, './param.in')
	except:
		print("Error: problem copying param.in from %s to run dir" % SCRIPT_PATH)
		raise

	try:
		_cmd_args = ["sed", "-i", "151s/ X / %i /" % n_days, "param.in"]
		subprocess.check_call(_cmd_args)
		if n_tracers > 0:
			_cmd_args = ["sed", "-i", "24s/ 0 / %i /" % n_tracers, "param.in"]
			subprocess.check_call(_cmd_args)
# TODO: Okay, this is stupid.  Fix this soon.
#			for i in range(n_tracers):
#				_line_number = 290
#				_cmd_args = ["sed", "-i", "%is/0/%i/" % (_line_number, i+1), 
#							 "param.in"]
#				subprocess.check_all(_cmd_args)
#				_line_number += 1
#				if i==2:
#					print("You must create tracer ouput files in param.in if"
#						  " you are using more than 2 tracers.\n")
#					break
	except:
		print("Error: problem editing param.in. Check file for errors.")
		raise

	try:
		shutil.copy("%s/../../bin/pelfe" % SCRIPT_PATH, "./pelfe")
	except:
		print("Error: Problem copying pelfe to run dir")
		raise

	try:
		shutil.copy('%s/selfe_run' % SCRIPT_PATH, './selfe_run')
	except:
		print("Error: Problem copying selfe_run to run dir")
		raise

	try:
		shutil.copy("%s/autocombine.pl", "./autocombine")
	except:
		print("Error: Problem copying selfe_run to run dir")
		raise

# If user wants an email message about run, do it.
# Maybe create an object for the selfe_run script 
	if email:
		try:
			_selfe_run = open("selfe_run~", "wb")
			_cmd_args = ["sed", "/#$ -S/{p;s/.*/#$ -M %s/;p;s/.*/#$ -m ae/;}"
						  % email, "selfe_run"]
			subprocess.check_call(_cmd_args, stdout = _selfe_run)
			_selfe_run.close()
			_cmd_args = ["mv", "-f", "selfe_run~", "selfe_run"]
			subprocess.check_call(_cmd_args)
		except:
			print("Error: problem editing selfe_run.  Check file for errors.")
			raise
# Name the run based on db and date
	try:
		_selfe_run = open("selfe_run~", "wb")
		_cmd_args = ["sed", "/#$ -S/{p;s/.*/#$ -N db%i_%i%i/;}" % (grid, 
				    startMonth, startYear), "selfe_run"]
		subprocess.check_call(_cmd_args, stdout = _selfe_run)
		_selfe_run.close()
		_cmd_args = ["mv", "-f", "selfe_run~", "selfe_run"]
		subprocess.check_call(_cmd_args)
	except:
		print("Error: problem editing selfe_run.  Check file for errors.")
		raise
# Clean up your mess
# subprocess.Popen(["rm", "-f", "*.new", "*.dat", "fort.*", "*.temp", \
# "date.in"], shell=True)
	os.system("rm -f *.new *.dat fort.* *.temp date.in")
	print ("Finished. Baseline run setup, specific implementation details not" 
           " handled by the script and must be made manually.\n")

if __name__ == '__main__':
# Create a parser for command line args and options
	usage = ("Usage: %prog [begin_date mm-dd-yyyy] [number of days in run] "
             "[grid [15(hi-res), 22, 26, 27(26 hi-res)]] [options]")

	parser = OptionParser(usage=usage)
	parser.add_option("-t", "--tracer", action="store", type="int",
					  dest="n_tracers", help="number of tracers")
	parser.add_option("-m", "--mail", action="store", type="string",
					  dest="email_addr", help="email to send run message")
	parser.set_defaults(n_tracers=0)
	(options, args) = parser.parse_args()
	n_tracers = options.n_tracers
# If email is not used, send it as empty string.  Hack to work until a more
# robust solution is created.
	if options.email_addr:
		email = options.email_addr
	else:
		email = ''

# grab start date, number of days, and db(grid) to usea
	if len(args) != 3:
		parser.error("Incorrect number of arguments")
	else:
		start = args[0]
		n_days = int(args[1])
		grid  = int(args[2])

#  If grid != 22 or 26 die - 27 is super secret code for high res version of 26
	if(grid!=15 and grid!=22 and grid!=26 and grid!=27):
		sys.exit("Grid options are: 15(hi-res), 22, 26, or 27(26 hi-res)")

# deal with start and ending days
	(month, day, year) = start.split('-')
	startDay   = int(day)
	startMonth = int(month)
	startYear  = int(year) 

# python style dates
	startDate = datetime.date(startYear, startMonth, startDay)

	selfeSetup(startDate, n_days, grid, n_tracers, email)
