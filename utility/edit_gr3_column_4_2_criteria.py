#!/usr/local/bin/python
""" This script is used to edit the fourth column in a *.gr3 file.

Edits gr3 column for between criteria_less and criteria_more
"""

import sys
import linecache

def change_gr3(old_file, new_file, criteria_value_less, criteria_value_more): 
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
		elems = int(_blah[0])
		nodes = int(_blah[1])
		for i in range(nodes):
			_blah = old_file.readline()
			_blah = _blah.split()
			depth = float(_blah[3])
			if (depth >= criteria_value_less) & (depth <= criteria_value_more):
				interp = depth/10 
			else:
				interp = depth 
			new_file.write("%s\t\t%s\t\t%s\t\t%.6e\n" % (_blah[0], _blah[1], _blah[2], interp))
		for i in range(elems):
			_blah = old_file.readline()
			new_file.write("%s" % _blah)
	except:
		raise
	finally:
		old_file.close()
		new_file.close()

if __name__ == '__main__':
	if len(sys.argv) != 5:
		sys.exit("Usage: %s [hgrid.gr3 path] [newFileName] [criteria value] " 
		         "[less than criteria] [greater than criteria]" % sys.argv[0])
	else:
		old_file = sys.argv[1]
		new_file = sys.argv[2]
		criteria_less = float(sys.argv[3])
		criteria_more= float(sys.argv[4])
	
	change_gr3(old_file, new_file, criteria_less, criteria_more) 
	
