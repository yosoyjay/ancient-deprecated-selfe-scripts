#!/usr/local/bin/python
""" This script is used to make all xyz values floats. 

"""

import sys
import linecache

def change_gr3(old_file, new_file): 
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
			new_file.write("%s\t\t%s\t\t%s\t\t%s\n" % (_blah[0], float(_blah[1]), float(_blah[2]), depth))
	except:
		raise
	finally:
		old_file.close()
		new_file.close()

if __name__ == '__main__':
	if len(sys.argv) != 3:
		sys.exit("Usage: %s [hgrid.gr3 path] [newFileName] " % sys.argv[0])
       	else:
		old_file = sys.argv[1]
		new_file = sys.argv[2]
	
	change_gr3(old_file, new_file)
	
