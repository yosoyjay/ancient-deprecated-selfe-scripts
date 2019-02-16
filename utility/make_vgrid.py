#!/usr/local/bin/python
""" Creates a s-level only vgrid.in file for use with SELFE 

    Only creates s-levels and assumes the max depth is < 100m.
    
    Example:
        ./make_vgrid.py 21 

        Creats a vgrid.in.new file with 21 s levels
"""
#-------------------------------------------------------------------------------
# Imports
#-------------------------------------------------------------------------------
import sys
import os

import numpy as np
#-------------------------------------------------------------------------------
# Constants 
#-------------------------------------------------------------------------------
ZSTD = [-5000, -2300, -1800, -1400, -1000, -770, -570, -470, -390,
		    -340,  -290,  -240,  -190,  -140,  -120,  -110, -105, -100]
NZLEVELS = 18

H_S = 100
H_C = 30.0
THETA_B = 1.0 
THETA_F = 10.0
NSLEVELS = 37

NLEVELS = 54 
#-------------------------------------------------------------------------------
# Classes and functions
#-------------------------------------------------------------------------------
class VGrid(object):
  """Vgrid related methods"""
  def __init__(self, vgridPath=None):
    """Creates vgrid object
  
    Args:
      vgridPath -- String of path to vgrid.in file to read in
    """
    self.zLevels = ZSTD 
    self.sLevels = []
    self.Levels = NLEVELS 
    self.nZLevels = NZLEVELS
    self.nSLevels = NSLEVELS 
    self.h_s = H_S 
    self.h_c = H_C 
    self.theta_b = THETA_B 
    self.theta_f = THETA_F 

    if vgridPath:
      self.read_Vgrid(vgridPath)

  def read_Vgrid(self, vgridPath):
    if not os.path.isfile(vgridPath): 
      raise Exception('File does not exist: '+vgridPath)
    f = open(vgridPath)
    # Header 
    self.nLevels, self.nZLevels, self.h_s = f.readline().split()
    # Z Levels
    f.readline()
    for line in xrange(zLevels+1):
      self.zLevels.append(f.readline())
    # S Levels
    f.readline()
    self.h_c, self.theta_b, self.theta_f = f.readline().split()
    for line in xrange(self.nLevels-self.nZLevels+1):
      self.sLevels.append(f.readline())
  
  def makeVGrid(self, fileName=None, zLevels=None, h_s=None, nSLevels=None, 
    h_c=None, theta_b=None, theta_f=None):
    """Creates a new vgrid.in file.
  
    Args:
      fileName -- String to path/filename of new vgrid.in file. 
            'vgrid.new' defaults
      zLevels -- Integer of the number of z-levels in new file
        18 default
      h_s -- Float of the transition depth from S to Z level
        100.0 default
      nSLevels -- Integer of the number of s-levels.
        37 default
      h_c -- Float of the depth for stretching
        30 default
      theta_b -- Float of the distribution of s-levels to bottom?
        0.7 default
      theta_f -- Float of distribution (depth of pycnocline?)
        10.0 default
    Returns:
      NOTHING
    Outputs:
      New vgrid file
    """
    if fileName != None:
      if not os.path.isfile(fileName):
        raise Exception('Unable to open file: '+fileName)
      file = fileName
    else:
      file = 'vgrid.in.new'
  
    if zLevels != None:
      self.zLevels = zLevels
    if h_s != None:
      self.h_s = h_s
    if nSLevels != None: 
      self.nSLevels = nSLevels
    if h_c != None:
      self.h_c = h_c
    if theta_b != None:
      self.theta_b = theta_b
    if theta_f != None:
      self.theta_f = theta_f      
 
    self.nLevels = self.nSLevels + self.nZLevels - 1
    self.sLevels = np.linspace(-1.0, 0.0, self.nSLevels)
  
    try:
      vgrid = open(file,'w')
    except:
      raise Exeption('Unable to write to file: '+file)

    vgrid.write('%d %d %.1f\n' % (self.nLevels, self.nZLevels, self.h_s))
    vgrid.write('Z levels\n')
    for l in xrange(self.nZLevels):
      vgrid.write('%d %d\n' % (l+1, self.zLevels[l]))
    vgrid.write("S levels\n")
    vgrid.write('%.2f %.2f %.2f\n' % (self.h_c, self.theta_b, self.theta_f))
    for l in xrange(self.nSLevels): 
      vgrid.write('%d %f\n' % (l+1, self.sLevels[l]))
    vgrid.close()
  
if __name__ == '__main__':
  if len(sys.argv) != 2:
    sys.exit("Usage: %s [number of S levels]" % sys.argv[0]); 
  else:
    vg = VGrid()
    vg.makeVGrid(nSLevels=int(sys.argv[1]))
