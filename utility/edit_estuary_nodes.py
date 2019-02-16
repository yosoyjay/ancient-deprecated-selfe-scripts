#!/usr/local/bin/python
""" This script is used to edit values that are either inside or outside the estuary

I made this to create windfactor.gr3 such that all the values in the estuary are given
value A and those outside the estuary are given value B.

Example:  I need nodes in estuary to be given value of 0 and nodes outside the estuary
are given a value 1.

./edit_estuary_nodes.py estuary.gr3 windfactor.gr3 0 1 
"""

import sys
import linecache

def edit_estuary_nodes(oldFile, newFile, outside, inside):
	try:
		oldFile = open(oldFile, "r")
		newFile = open(newFile, "w")
	except:
		raise

	try:
# Copy the header
		_tmp = oldFile.readline()
		newFile.write("%s" % _tmp)
		_tmp = oldFile.readline()
		newFile.write("%s" % _tmp)
# Get number of nodes.  
# If in estuary [_tmp(3) == 1] give in value
# else give out value
		_tmp = _tmp.split()
		nElems = int(_tmp[0])
		nNodes = int(_tmp[1])
		for i in range(nNodes):
			_tmp = oldFile.readline()
			_tmp = _tmp.split()
			locFlag = int(float(_tmp[3]))
			if locFlag == 1:
				newFile.write("%s\t\t%s\t\t%s\t\t%s\n" % (_tmp[0], _tmp[1], _tmp[2], inside))
			else:
				newFile.write("%s\t\t%s\t\t%s\t\t%s\n" % (_tmp[0], _tmp[1], _tmp[2], outside))
		for i in range(nElems):
			_tmp = oldFile.readline()
			newFile.write("%s" % _tmp)
	except:
		raise
	finally:
		oldFile.close()
		newFile.close()

if __name__ == '__main__':
	if len(sys.argv) != 5:
		sys.exit("Usage: %s [estuary.gr3 path] [newFileName] [outsideEstuary] [insideEstuary]"
		         % sys.argv[0])
	else:
		oldFile = sys.argv[1]
		newFile = sys.argv[2]
		outside = sys.argv[3]
		inside  = sys.argv[4]

	edit_estuary_nodes(oldFile, newFile, outside, inside)	
