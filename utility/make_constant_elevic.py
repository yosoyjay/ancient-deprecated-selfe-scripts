#!/usr/local/bin/python
""" This script is used to edit the fourth column in a *.gr3 file.

This script was created to create an elev.ic file based on an hgrid.gr3.
The water column will have a constant depth of 20m.

Example: I need a constant depth of 20m everywhere.

Input: ./thisScript.py 20

Output: elev.ic
"""

import sys
import linecache

def change_gr3(old_file, new_file, constDepth): 
	try:
		old_file = open(old_file, 'r')
	except:
		raise

	try:
		new_file = open(new_file, 'w')
	except:
		raise

	try:
# Copy the header from the old file
	   	_blah = old_file.readline()
		new_file.write("%s" % _blah)
		_blah = old_file.readline()
		new_file.write("%s" % _blah)
# Get number of nodes in the file and iterate over that number
# Get the 4th column value and determine what the new value should be
# and write it back to the new file.
		_blah = _blah.split()
		nodes = int(_blah[1])
		for i in range(nodes):
			_blah = old_file.readline()
			_blah = _blah.split()
			depth = float(_blah[3])
			depth = constDepth - depth
			new_file.write("%s\t%s\t%s\t%f\n" % (_blah[0], _blah[1], _blah[2], depth))
	except:
		raise
	finally:
		old_file.close()
		new_file.close()

if __name__ == '__main__':
	if len(sys.argv) != 4:
		sys.exit("Usage: %s [hgrid.gr3 path] [newFileName] [imposed depth] " 
		          % sys.argv[0])
	else:
		old_file = sys.argv[1]
		new_file = sys.argv[2]
		constDepth = float(sys.argv[3])
	
	change_gr3(old_file, new_file, constDepth) 
	
