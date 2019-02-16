#!/usr/local/bin/python
""" Creates a nudging file for an open channel hgrid  

    Example:
        ./make_nudge.py hgrid.gr3 20000 

        Creates a nudge.in file where the 4th column linearly diminishes from
        1000 to 0 over the first 20000 units of the grid.
"""
import sys
from numpy import array, empty

def make_nudge(gridPath, endXPos):
  # Constants
  highVal = 100
  lowVal = 0

  # Load hgrid and get nodes
  [nodes, elements] = _loadGrid(gridPath)
  nudge = open("nudge.in", "w")

  # Calculate delta x grid --- Assummed to be constant
  x0 = nodes[0][0]
  idx = 1
  while x0 == nodes[idx][0]:
      idx = idx + 1
  x1 = nodes[idx][0]
  dx = x1-x0

  # Calculate slope
  m = (highVal - lowVal)/endXPos

  # Write nudging file 
  nudge.write("Nudge from %f to %f\n" % (highVal, lowVal))
  nudge.write("%d %d\n" % (elements.shape[0], nodes.shape[0]))
  for node in range(nodes.shape[0]):
    if nodes[node][0] < endXPos:
        nudge.write("%d %f %f %f\n" % (node+1,nodes[node][0],nodes[node][1],highVal-dx*m*(nodes[node][0]/dx))) 
    else:
        nudge.write("%d %f %f %f\n" % (node+1,nodes[node][0],nodes[node][1],0))
 
def _loadGrid(gr3Path):
    """ Loads a *.gr3 file into a numpy array for nodes and elements """
    try:
        gr3 = open(gr3Path, "r")
        try:
            # read header, elements, and nodes
            gr3.readline()
            _input = (gr3.readline()).split()
            nNodes = int(_input[1])
            nElems = int(_input[0])

            # read in all nodes
            nodes = empty((nNodes,3))
            for i in range(nNodes):
                _input = (gr3.readline()).split()
                nodes[i] = [float(_input[1]),float(_input[2]),float(_input[3])]

            # read in elements
            elements = empty((nElems,3))
            for i in range(nElems):
                _input = (gr3.readline()).split()
                elements[i] = [int(_input[2])-1,int(_input[3])-1,int(_input[4])-1]
        finally:
            gr3.close()
    except IOError:
        print("Error: Problem loading *.gr3 file check path")
        raise

    return nodes, elements
                


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit("Usage: %s [grid] [end]" % sys.argv[0]); 
    else:
        grid = sys.argv[1] 
        end = float(sys.argv[2])

    make_nudge(grid, end)

