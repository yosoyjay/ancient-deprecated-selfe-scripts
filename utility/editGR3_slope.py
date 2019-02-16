""" This script is used to edit the fourth column in a *.gr3 file.

Use this script to apply a linear gradient along the x-direction.
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
START = 30000
END = 80000
MINVAL = 30 
MAXVAL = 0 

#-------------------------------------------------------------------------------
# Functions 
#-------------------------------------------------------------------------------
def calcDepth(x, start, end, minVal, maxVal):
  """Calculates linear slope and applies gradient based on x. 

  No min or max limits applied.
  
  Args:
    x - Float of x positions (or value) 
    start - Float of where to start applying linear gradient
    end - Float of where to stop applying linear gradient
    minVal - Float of value in gradient at start point
    maxVal - Float of value in gradient at end point 
  """
  slope = float(maxVal - minVal)/float(end - start)
  if x < start:
    val = minVal
  elif x > end:
    val = maxVal
  else:
    val = slope*(x - start) + minVal
  return val

def applyFunctionToGrid(hgridFile, newFile):
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
  nodes = int(header.rstrip().split()[1])
  fNew.write(header)

  # Get estuary flags, read in nodes, calc value, and write out to new file
  for i in range(nodes):
    line = fOld.readline().rstrip().split()
    x = float(line[1])
    y = float(line[2])
    depth = float(line[3])
    value = calcDepth(x, START, END, MINVAL, MAXVAL)
    fNew.write('%d\t\t%.2f\t\t%.2f\t\t%.6f\n' % (i+1, x, y, value))

  # Copy element table
  for i in xrange(elements):
    line = fOld.readline()
    fNew.write(line)

  # Copy open boundaries
  tmp = fOld.readline()
  fNew.write(tmp)
  nOpenBnds = int(tmp.rstrip().split()[0])
  fNew.write(fOld.readline())
  for openBnd in range(nOpenBnds):
    tmp = fOld.readline()
    fNew.write(tmp)
    nOpenBndNodes = int(tmp.rstrip().split()[0])
    for i in range(nOpenBndNodes):
      fNew.write(fOld.readline()) 
  # Copy closed boundaries
  tmp = fOld.readline()
  fNew.write(tmp)
  nOpenBnds = int(tmp.rstrip().split()[0])
  fNew.write(fOld.readline())
  for openBnd in range(nOpenBnds):
    tmp = fOld.readline()
    fNew.write(tmp)
    nOpenBndNodes = int(tmp.rstrip().split()[0])
    for i in range(nOpenBndNodes):
      fNew.write(fOld.readline()) 
  # Copy island boundaries if they exist
  try:
    tmp = fOld.readline()
    fNew.write(tmp)
    nOpenBnds = int(tmp.rstrip().split()[0])
    fNew.write(fOld.readline())
    for openBnd in range(nOpenBnds):
      tmp = fOld.readline()
      fNew.write(tmp)
      nOpenBndNodes = int(tmp.rstrip().split()[0])
      for i in range(nOpenBndNodes):
        fNew.write(fOld.readline()) 
    fOld.close()
    fNew.close()
  except:
    fOld.close()
    fNew.close()


if __name__ == '__main__':
  if len(sys.argv) != 3:
    sys.exit("Usage: %s [hgrid.gr3 path] [newFileName]" % sys.argv[0]) 
  else:
    hgridFile = sys.argv[1]
    newFile = sys.argv[2]
       
  applyFunctionToGrid(hgridFile, newFile)
