#!/usr/local/bin/python
"""
This script generates the tidal amplitudes and phases for bctides
based upon Joe Cho's guide.  

It's more or less hardcoded, kludgy, and crappy.
"""
import sys
import os
import datetime
import shutil
import subprocess

def genTidalData(startDate, grid):
# Some constants like directories and such
	HOME_PATH    = "/home/workspace/project/lopezj/"
	SCRIPT_PATH  = "/home/workspace/project/lopezj/scripts/selfe_setup/"
	BIN_PATH     = "/home/workspace/project/lopezj/bin/selfe_setup/"
	TIDES_PATH   = "/home/workspace/project/lopezj/scripts/selfe_setup/tides/"

# Create new files to get tidal constituents and strip away 
# what is not needed.
# Get number of elements and nodes from hgrid.gr3 files to be used
# TODO: Fix this crap - Make a grid object
	header = []
	nodes  = []
	table  = []
	boundary = []
	try:
		hgrid = open("hgrid.gr3", "r")
		try:
# Get header info, nNodes, and nElements
			header.append(hgrid.readline())
			_input = hgrid.readline()
			header.append(_input)
			_input = _input.split(" ")
			nNodes = int(_input[1])
			nElements = int(_input[0])
# Get all node info and element info
			for i in range(nNodes):
				nodes.append(hgrid.readline())
			for i in range(nElements):
				table.append(hgrid.readline())
# Get boundary info by reading rest of the file and save in a list and file
			while True:
				_input = hgrid.readline()
				if len(_input) == 0:
					break
				boundary.append(_input)
		finally:
			hgrid.close()
	except IOError:
		print("Error: problem reading in hgrid.gr3")
		raise	

# Write the files that are used later.  This crap needs to be fixed.
# Write fort.14
	try:
		fort14 = open("fort.14", "wb")
		try:
			for i in range(nNodes):
				fort14.write("%s" % nodes[i])
		finally:
			fort14.close()
	except IOError:
		print("Error: problem writing fort.14")
		raise

# Write boundary.new
	try:
		bndry = open("boundary.new", "wb")
		try:
			for i in range(len(boundary)):
				bndry.write(boundary[i])
		finally:
			bndry.close()
	except IOError:
		print("Error: problem writing boundary.new")
		raise

# Stupidness to read in boundaries separated.  Please fix this soon.
# Will be fixed with creation of grid objects
	try:
		bndry = open('boundary.new', 'rb')
		try:
			_input = bndry.readline()
			bndry.readline()
			_input = _input.split(" ")
			nBndry = int(_input[0])
			bndryList = []
			for eachBndry in range(nBndry):
				_input = bndry.readline()
				_input = _input.split(" ")
				nBndryNodes = int(_input[0])
				nodeList = []
				for eachNode in range(nBndryNodes):
					_input = int(bndry.readline())
					nodeList.append(_input)
				bndryList.append(nodeList)
		finally:
			bndry.close()
	except IOError:
		print("Error: Problem reading in boundary.new")
		raise

# Run spcs.pl to get fort.14.ll
	sys.path.append("%sbin/pexpect/" % HOME_PATH)
	import pexpect
	child = pexpect.spawn("%sspcs.pl" % TIDES_PATH, timeout=90)
	child.expect("Enter input file name to be converted: ")
	child.sendline("fort.14")
	child.expect("Enter output file name: ")
	child.sendline("fort.14.ll")
	child.expect("\n")
	child.expect("\n")
	child.expect("\n")
	child.expect("\n")
	child.expect("Enter data type: ")
	child.sendline("2")
	child.expect("\n")
	child.expect("\n")
	child.expect("\n")
	child.expect("\n")
	child.expect("Enter data type: ")
	child.sendline("1")
#	child.expect("Enter a digit [0-8] for the precision 'proj' uses for output: ")
#   This line above hangs... I thing the ' need to be escaped.
	child.sendline("8")
	child.expect("\n")
	child.expect("\n")
	child.expect("\n")
	child.expect("\n")
	child.expect("Enter style number: ")
	child.sendline("1")
	child.expect("Enter state plane coordinate system zone number: ")
	child.sendline("3601")
	child.expect("Enter NADCON correction region: ")
	child.sendline("wo")
	child.expect("\n")
	child.expect("\n")
	child.expect("\n")
	child.expect("Enter units type: ")
	child.sendline("2")
	child.expect("Done")

# Fix up new fort.14.ll file with header, element table, and boundary info
# to retain same structure as *.gr3 files.
	fort14ll = []
	try:
		fort14 = open("fort.14.ll", "rb")
		try:
			while True:
				_input = fort14.readline()
				if len(_input) == 0:
					break
				fort14ll.append(_input)
		finally:
			fort14.close()
	except IOError:
		print("Error reading in fort.14.ll")
		raise

	try:
		fort14 = open("fort.14.ll", "wb")
		try:
			for i in range(len(header)):
				fort14.write("%s" % header[i])
			for i in range(len(fort14ll)):
				fort14.write(fort14ll[i])
			for i in range(len(table)):
				fort14.write(table[i])
			for i in range(len(boundary)):
				fort14.write(boundary[i])
		finally:
			fort14.close()
	except IOError:
		print("Error writing to fort.14.ll")
		raise

# Run ecp to get *.nos8 file and add boundary info
	try:
		_cmd_args = "%secp" % BIN_PATH
		subprocess.check_call(_cmd_args)
	except IOError:
		print("Problem creating fort.14.no8")
		raise
	
# needed for intel_get to get ap.dat files
	os.system("ln -fs %sintel_deg    ./" % TIDES_PATH)
	os.system("ln -fs %sedpac2xy.gr3 ./" % TIDES_PATH)
	os.system("ln -fs %steanl.tct    ./" % TIDES_PATH)
	os.system("cp fort.14.nos8 nos8")				# hack to get it working.

# Hard coded to just get the first two boundaries - Strait of Georgia and Ocean
# TODO: Fix this crap, it's really, really bad
# It's even worse now that I put in an if statement for db15!
# Make grid objects
	fakeIter = 0
	for eachBndry in bndryList:
		if grid == 15:
			if fakeIter == 1:
				break
		if fakeIter == 2:
			break
		destination = open('nos8_%i.new' % fakeIter, 'wb')
		shutil.copyfileobj(open('nos8', 'rb'), destination)
		destination.write('%i\n' % len(bndryList))  		   	# of boundaries
		destination.write('%i\n' % len(bndryList))  		   	# doesn't matter 
		destination.write('%i\n' % len(bndryList[fakeIter]))  	# nodes in bndry
		for eachNode in eachBndry:
			destination.write('%i\n' % eachNode)			   	# write node number
		destination.close()

		try:
			_cmd_args = ["cp", "-f", "nos8_%i.new" % fakeIter, "fort.14.nos8"]
			subprocess.check_call(_cmd_args)
		except :
			print("Error: Problem copying nos8_i% to fort.14.nos8" % fakeIter)
			raise

# Run genbc to get *.sta 
		try:
			_cmd_args = ["%sgenbcs" % TIDES_PATH]
			subprocess.check_call(_cmd_args)
		except :
			print("Error: Problem with call to genbcs.f")
			raise
	
# run intel_deg to get ap.dat
		child2 = pexpect.spawn("./intel_deg", timeout=90)
		child2.expect("Grid file : ")
		child2.sendline("edpac2xy.gr3")
		child2.expect("Name of file with location of stations: ")
		child2.sendline("fort.14.sta")
		child2.expect("TEANL file : ")
		child2.sendline("teanl.tct")
		for i in range(len(bndryList[fakeIter])):
			child2.expect("\n")
		child2.expect("Output file: ")
# I want the files to start at 1 not zero
		fakeIter = fakeIter+1				
		child2.sendline("ap_%i.dat" % fakeIter)

# Copy tide_colu.com, change date, and run tide1_e.
	try:
		_cmd_args = ["cp", "%s/tide_colu.com" % TIDES_PATH, "./"]
		subprocess.check_call(_cmd_args)
		_cmd_args = ["cp", "%s/fort.8" % TIDES_PATH, "./"]
		subprocess.check_call(_cmd_args)
	except :
		print("Error: Problem creating symlinks required for tide1_e\n"
			  "Check for tide_colu.com and fort.8")
		raise

	try:
		year = startDate[0] - 2000	# Used for weird date required for tide1_e
		_cmd_args = ["sed", "-i", "1s/dateHolder/08%s%s%i20/" %
			(str(startDate[2]).zfill(2), str(startDate[1]).zfill(2), year), 
			"tide_colu.com"]
		subprocess.check_call(_cmd_args)
	except :
		print("Error: Problem editing tide_colu.com, check file to ensure "
		      "correct parameters.")
		raise

	os.system("%s/tide1_e < tide_colu.com > tide_colu.out" % BIN_PATH)
#	try:
#		_file_in = open("tide_colu.com", "rb")
#		_file_out = open("tide_colu.out", "wb")
#		try:
#			_cmd_args = ["%s/tide1_e" % BIN_PATH] 
#			subprocess.check_call(_cmd_args, stdin=_file_in, stdout=_file_out)
#		except :
#			print("Error: Problem calling tide1_e.f")
#			raise
#		finally:
#			_file_in.close()
#			_file_out.close()
#	except IOError:
#		print("Error: Problem opening tide_colu.com or tide_colu.out for "
#		      "tide1_e.f")
#		raise

# Clean up temp files
	try:
		_cmd_args = ["rm", "-f", "nos8", "intel_deg", "ssh.in"]
		subprocess.check_call(_cmd_args)
	except :
		print("Error: Problem deleting temporary files")
		raise

if __name__ == '__main__':
# grab start and number of days 
	if len(sys.argv) != 3:
		sys.exit("usage: %s [begin_date mm-dd-yyyy] [grid [15,22,26,27]" 
		         % sys.argv[0])
	else:
		start = sys.argv[1]
		grid = int(sys.argv[2])

# deal with start days
	(month, day, year) = start.split('-')
	startDay   = int(day)
	startMonth = int(month)
	startYear  = int(year) 

# Python style dates for script.
	startDate = datetime.date(startYear, startMonth, startDay)
	startDate = startDate.timetuple()

	genTidalData(startDate, grid)
