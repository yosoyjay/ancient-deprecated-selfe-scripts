#!/usr/local/bin/python
""" Creates a *.bp transect file in a straight line for open channel 

    
    Example:
        ./make_bp.py 0 40000 200 25 

        Creats a new.bp file with build points every 200 meters from 0 to 40000
        with a constant y value of 25.
"""
import sys
from math import ceil

def make_bp(start, end, spacing, y):
    try:
        bp = open('new.bp','w')
        bp.write("Open channel transect from %d m to %d m every %d m "
                 "at %d meters in y\n" % (start,end,spacing,y))
	# Number of build points
	nBps = int(ceil((end-start)/spacing))
        bp.write("%d\n" % nBps)
        for i in range(1, nBps+1):
            bp.write("%d\t%f\t%f\t1.0\n" % (i, start+(spacing*(i-1)), y))       
    except:
        raise
    finally:
        bp.close()

if __name__ == '__main__':
    if len(sys.argv) != 5:
        sys.exit("Usage: %s [start] [end] [spacing] [y]" % sys.argv[0]); 
    else:
        start = float(sys.argv[1])
        end = float(sys.argv[2])
        spacing = float(sys.argv[3])
        y = float(sys.argv[4])

    make_bp(start, end, spacing, y)

