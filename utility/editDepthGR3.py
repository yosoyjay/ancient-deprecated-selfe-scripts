#!/usr/local/bin/python
""" This script is used to edit the fourth column in a *.gr3 file.

I'm using this to change an hgrid.gr3 to lqk.gr3 for interpolation in the 
ELM backtracking.  It changes the values in the fourth column based on 
a less than input or greater than input criteria.

value_less and value_more written as integers. 

if $4 <= criteria_value, $4 = value_less
else $4 = value_more

Example:  I need nodes with depth <= 100 to set $4 to 2 and nodes with
depth > 100 to set $4 to 1.
criteria_value = 100
value_less = 2
value_less = 1
"""

import sys
import linecache

def change_gr3(old_file, new_file, criteria_value, value_less, value_more): 
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
                        if depth <= criteria_value:
                                interp = value_less
                        else:
                                interp = value_more
#                                interp = depth
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
        if len(sys.argv) != 6:
                sys.exit("Usage: %s [hgrid.gr3 path] [newFileName] [criteria value] " 
                         "[less than criteria] [greater than criteria]" % sys.argv[0])
        else:
                old_file = sys.argv[1]
                new_file = sys.argv[2]
                criteria = float(sys.argv[3])
                value_less = float(sys.argv[4])
                value_more = float(sys.argv[5])
        
        change_gr3(old_file, new_file, criteria, value_less, value_more)
        
