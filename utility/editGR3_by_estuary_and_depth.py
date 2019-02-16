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

#-------------------------------------------------------------------------------
# Imports 
#-------------------------------------------------------------------------------
import sys
import linecache
import numpy as np

#-------------------------------------------------------------------------------
# Constants 
#-------------------------------------------------------------------------------
# Lots of numbers saved here to remember how I derived it
CRITICAL_D = 5 
MIN_VAL = 0.0005
MAX_VAL = 0.005
# Just to remind myself in case I need it at some point. These are eyeballed
# values.
MOUTH_X = 333785 
ASTORIA_X = 349334 
TONGUE_POINT_X = 358310
B_SLOPE = (MAX_VAL - MIN_VAL) / (TONGUE_POINT_X - ASTORIA_X) 
A_SLOPE = (MAX_VAL - MIN_VAL) / (TONGUE_POINT_X - MOUTH_X) 

#-------------------------------------------------------------------------------
# Functions 
#-------------------------------------------------------------------------------
def calcBelow(x):
  """Calculates and returns a value based on x value (or whatever) """
  return min(max(B_SLOPE*(x - ASTORIA_X) + MIN_VAL, MIN_VAL), MAX_VAL)

def calcAbove(x):
  """Calculates and retursn a value based on x value (or whatever) """
  return min(max(A_SLOPE*(x - MOUTH_X) + MIN_VAL, MIN_VAL), MAX_VAL) 

def getEstuary(estuaryFile):
  """Gets estaury value and returns array with estuary values [0:no, 1:yes]

  Args:
    estuaryFile -- Path to estuary.gr3 file

  Returns:
    estuaryFlags -- Numpy array of flags indicating whether node is in estuary 
                    or not.
  """
  try:
    f = open(estuaryFile, 'r')
  except IOError:
    raise Exception( 'File does not exist:'+estuaryFile )

  header = f.readline()
  header = f.readline()
  elems, nodes = header.split()   
  elems = int(elems)
  nodes = int(nodes)

  estuary = np.zeros([nodes, 1])
  for i in range(nodes):
    tmp = f.readline().split()
    if int(float(tmp[3])) != 0:
      estuary[i] = 1
    else:
      estuary[i] = 0

  f.close()
  return estuary

def applyFunctionToGrid(hgridFile, estuaryFile, newFile):
  """Gets hgrid and creates a newFile based on a function.

  Args:
    hgridFile -- String path to hgrid.gr3
    estuaryFile -- String of path to estuary.gr3
    newFile -- String of path to new file to create based on function
  Returns:
    Nothing
  """
  try:
    fOld = open(hgridFile, 'r')
  except IOError:
    raise Exception( 'File does not exist:'+hgridFile )

  try:
    fNew = open(newFile, 'w')
  except IOError:
    raise Exception( 'Cannot open file:'+newFile )

  # Copy the header
  header = fOld.readline()
  fNew.write(header)
  header = fOld.readline()
  elements = int(header.rstrip().split()[0])
  fNew.write(header)

  # Get estuary flags, read in nodes, calc value, and write out to new file
  estuary = getEstuary(estuaryFile)
  for i in range(estuary.shape[0]):
    line = fOld.readline().rstrip().split()
    x = float(line[1])
    y = float(line[2])
    depth = float(line[3])
    if estuary[i] == 1:
      if depth > CRITICAL_D:
        value = calcBelow(x)
      else:
        value = calcAbove(x)
    else:
      value = MIN_VAL
    fNew.write('%d\t\t%.2f\t\t%.2f\t\t%.6f\n' % (i+1, x, y, value))

  # Copy element table
  for i in xrange(elements):
    line = fOld.readline()
    fNew.write(line)

  # Do not copy boundaries, islands, etc.
  fOld.close()
  fNew.close()


if __name__ == '__main__':
  if len(sys.argv) != 4:
    sys.exit("Usage: %s [hgrid.gr3 path] [estuary.gre path] [newFileName]" % sys.argv[0]) 
  else:
    hgridFile = sys.argv[1]
    estuaryFile =  sys.argv[2]
    newFile = sys.argv[3]
       
  applyFunctionToGrid(hgridFile, estuaryFile, newFile)
