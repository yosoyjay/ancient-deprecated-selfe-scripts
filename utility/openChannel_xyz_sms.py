#!/usr/local/bin/python
""" This script is used to create a rectangular prism grid for SMS. 

Input:
1. x - Length along x-axis in meters
2. y - Length along y-axis in meters
3. z - Length along z-axis in meters
4. x_div - Number of segments along x-axis
5. y_div - Number of segments along y-axis
6. slope - Delta z / Delta x

jlopez - 9/1/2011 
"""

import math
import sys

def create_xyzFile(x, y, z, x_div, y_div, slope):

# Create new file. Name based on dimensions.
	try:
		_file = "%d_%d_%d_%d_%d.xyz" % (x, y, z, x_div, y_div)
		xyz_file = open(_file, 'w')
	except:
		raise

	try:
# Write header and loop over the x,y dimensions. Slope can vary by x-coordinate. 
		xyz_file.write("x:%d y:%d z:%d x-segments:%d y-segments:%d "
			           "\n" % (x, y, z, x_div, y_div))

		x_skip = float(x)/x_div
		y_skip = float(y)/y_div

		x_itr = int(math.ceil(x/x_skip))
		y_itr = int(math.ceil(y/y_skip))

		for i in range(x_div+1):
			for j in range (y_div+1):
				xyz_file.write("%f %f %f\n" % (i*x_skip,j*y_skip,z+(slope*i*x_skip)))
	
	except:
		raise
	finally:
		xyz_file.close()

if __name__ == '__main__':
	if len(sys.argv) != 7:
		sys.exit("Usage: %s [len of x dim] [len of y dim] [depth of z dim]"
				 " [# x segs] [# y segs] [slope]" % sys.argv[0])
	else:
		x = int(sys.argv[1])
		y = int(sys.argv[2])
		z = int(sys.argv[3])
		x_div = int(sys.argv[4])
		y_div = int(sys.argv[5])
		slope = float(sys.argv[6])

	create_xyzFile(x, y, z, x_div, y_div, slope)

