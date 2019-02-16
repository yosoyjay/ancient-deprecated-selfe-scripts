#!/usr/local/bin/python
""" Creates bctides.in for a selfe run

This file must be run in the sequence set forth in selfeSetup.py due 
this script extracting data from files generated in other scripts.
Crappy design, I know. 

Designed to be used with selfe V3.1d
"""
import os
import sys
import datetime
import linecache

def genBCTides( startDate, grid ):

# Constitutents used for NCOM forcing - Totally hard coded
	CONSTITS  = ['Z0','O1', 'K1', 'Q2', 'P1', 'K2', 'N2', 'M2', 'S2']
	LINE_NUMS = [ 1,   13,   21,   11,   19,   43,   32,   36,   41 ]  
	N_CONSTITS = 9

# Boundary conditions - Totally hard coded
# Look up selfe user manual for specifics of what the value mean
	if grid == 15:
		BNDRY_NAMES = [ 'Beaver', 'Ocean' ]
		bndryConds = {}
		bndryVals  = []

		_conds = '65 3 0 0 0 !Ocean'
		bndryConds[ BNDRY_NAMES[0] ] = _conds
		bndryVals.append([])

		_conds = '4 0 1 1 2 !Beaver Army'
		bndryConds[ BNDRY_NAMES[1] ] = _conds
		bndryVals.append( ['0.7\n', '0.0\n', '0.1\n'] )

		nNCOMBndry = 1

	elif grid == 22:
		BNDRY_NAMES = [ 'Georgia', 'Ocean', 'Beaver', 'Fraser' ]
		bndryConds = {}
		bndryVals  = []

		_conds = '3 3 0 0 0 !Georgia' 
		bndryConds[ BNDRY_NAMES[0] ] = _conds
		bndryVals.append([])

		_conds = '88 3 0 0 0 !Ocean'
		bndryConds[ BNDRY_NAMES[1] ] = _conds
		bndryVals.append([])

		_conds = '5 0 1 1 2 !Beaver Army'
		bndryConds[ BNDRY_NAMES[2] ] = _conds
		bndryVals.append( ['0.7\n', '0.0\n', '0.1\n'] )

		_conds = '3 0 2 3 3 !Fraser '
		bndryConds[ BNDRY_NAMES[3] ] = _conds
		bndryVals.append( ['0.0\n', '0.0\n', '0.0\n'] )

		nNCOMBndry = 2

	elif grid == 26:
		BNDRY_NAMES = ['Georgia', 'Ocean', 'Bonneville', 'Williamette','Fraser']
		bndryConds = {}
		bndryVals  = []

# Note: Georgia and the Ocean are switched in db26
		_conds = '3 3 0 0 0 !Georgia' 
		bndryConds[ BNDRY_NAMES[1] ] = _conds
		bndryVals.append([])

		_conds = '88 3 0 0 0 !Ocean'
		bndryConds[ BNDRY_NAMES[0] ] = _conds
		bndryVals.append([])

		_conds = '6 0 1 1 2 !Bonneville'
		bndryConds[ BNDRY_NAMES[2] ] = _conds
		bndryVals.append( ['0.7\n', '0.0\n', '0.1\n'] )

		_conds = '3 0 1 1 2 !Willamette'
		bndryConds[ BNDRY_NAMES[3] ] = _conds
		bndryVals.append( ['0.7\n', '0.0\n', '0.1\n'] )

		_conds = '3 0 2 3 3 !Fraser'
		bndryConds[ BNDRY_NAMES[4] ] = _conds
		bndryVals.append( ['0.0\n', '0.0\n', '0.0\n'] )

		nNCOMBndry = 2
	else:
		print 'Invalid grid type selected'	

# Get information from boundary.new about boundaries that use ncom
# nNodesBndry - Number of nodes at a boundary
# openBndry   - Number of open boundaries
	try:
		boundary = open('boundary.new', 'rb')
	except IOError:
		print("Error: problem reading from boundary.new")
		raise
	else:
		try:
			nNodesBndry = []								
			_input = boundary.readline() 
			_input = str(_input).split()
			nOpenBndry = int(_input[0])							
			boundary.readline()

			for nBndry in range(nNCOMBndry):
				_input = boundary.readline()
				_input = str(_input).split()
				nNodesBndry.append(int(_input[0]))
				for i in range(nBndry):	
					boundary.readline()
		finally:
			boundary.close()

# Get information from ap_?.out 
# For each constituent of each boundary there is:
# nNodes = number of nodes on boundary only on NCOM effected boundaries
# fFreq  = forcing frequency of each contituent 
# nodeAP = list of nodal amplitude and phase for each node along a boundary 
	boundary = []
	for i in range(1, nNCOMBndry+1):
		_thisBndry = []
		try:
			apFile = open('ap_%i.dat' % i, 'rb')
		except IOError:
			print "Error: unable to open ap_%i.dat" % i
			raise
		else:
			try:
				apFile.readline()								# Dump header
				nNodes = apFile.readline()
				nNodes = int(nNodes.rstrip('\n'))				# # of nodes 
				apFile.readline()								# Dump # of constits
				for eachConstituent in range(N_CONSTITS):
					_constit = {}								# Amp, phase, nNodes 
					_constit['nNodes'] = nNodes					# Number of nodes
					_input = apFile.readline()  				# Forcing frequency 
					_input = _input.rstrip('\n')				
					_constit['fFreq']  = float(_input)  	    
					apFile.readline()							# Dump constit name 
					_nodeData = []
					for allNodes in range( nNodes ):
						_nodeData.append(apFile.readline())		# Append node to list
					_constit['nodeAP'] = _nodeData				# Add amp and phase 
					_thisBndry.append(_constit)					# Append constituent 
				boundary.append(_thisBndry)						# Append bndary to list
			finally:
				apFile.close()

# Read in data from Z0.out
	try:
		z0File = open("Z0.out", "rb")
	except IOError:
		print("Error: unable to open Z0.out")
		raise
	else:
		try:
			z0 = []
			for i in range(nNCOMBndry):
				z0File.readline()								# Throw out header
				_bndry = []
				for i in range(boundary[i][0]['nNodes']):		# Get Z0 for all nodes 
					_bndry.append(z0File.readline())			# Just copy it 
				z0.append(_bndry)
		finally:
			z0File.close()

# read in data from tide_colu.out - totally hardcoded line numbers
# to be read in tide_colu.out the output from intel_deg
	constituents = []
	for i in range(N_CONSTITS):
		_constit = [CONSTITS[i]]						# Get name
		_constit.append(boundary[0][i]['fFreq'])  	# Copy forcing frequency
		_input = linecache.getline('tide_colu.out', LINE_NUMS[i])  
		_input = _input.split()
		_constit.append(float(_input[1]))			# Nodal factor
		_constit.append(float(_input[2]))			# Earth equilibrium
		constituents.append(_constit)				# Add it to the list

# Time to write the file
	try:
		bctides = open('bctides.in', 'wb')
		bctides.write('%s/%s/%i %s:%s:%s PST\n' % (str(startDate[1]).zfill(2),     
			str(startDate[2]).zfill(2), startDate[0], str(startDate[3]).zfill(2),  
			str(startDate[4]).zfill(2), str(startDate[5]).zfill(2)))											   
		bctides.write('%i %i\n' % (0, 40))  		# some crap  
		bctides.write('%i\n' % N_CONSTITS) 			# of constituients
		for i in range(N_CONSTITS):
			bctides.write('%s\n' % CONSTITS[i])
			bctides.write('%f %f %f\n' % (constituents[i][1], constituents[i][2], 
						  constituents[i][3]))
		bctides.write('%i\n' % nOpenBndry)			# of open boundaries

# Loop over all the boundaries and write the values
		for i in range(nOpenBndry):
			if i < nNCOMBndry:
				bctides.write('%s\n' % bndryConds[BNDRY_NAMES[i]])	# Bndry conds
				bctides.write('%s\n' % CONSTITS[0])					# Z0
				for j in range( boundary[i][0]['nNodes'] ): 		# Write Z0 
					bctides.write('%s' % z0[i][j] )
				for j in range(1, N_CONSTITS):						# Write others
					bctides.write('%s\n' % CONSTITS[j])
					for k in range( boundary[i][j]['nNodes'] ): 	# over all nodes 
						bctides.write('%s' % boundary[i][j]['nodeAP'][k] )
			else:
				bctides.write('%s\n' % bndryConds[BNDRY_NAMES[i]])	# Bndry conds 
				for j in range( len(bndryVals[i]) ):				# Constants 
					bctides.write('%s' % bndryVals[i][j])
	except IOError:
		print("Error: Problems writing to bctides.in")
		raise
	finally:
		bctides.close()
	
	os.system("rm boundary.new Z0.out tide_*")

if __name__ == '__main__':
# grab start and number of days 
	if len(sys.argv) != 3:
		sys.exit("usage: %s [begin_date mm-dd-yyyy] [grid 15, 22, or 26]" %sys.argv[0])
	else:
		start = sys.argv[1]
		grid  = int(sys.argv[2])

# deal with start days
	(month, day, year) = start.split('-')
	startDay   = int(day)
	startMonth = int(month)
	startYear  = int(year) 

# Python style dates for script.
	startDate = datetime.date(startYear, startMonth, startDay)
	startDate = startDate.timetuple()

	genBCTides(startDate, grid)
