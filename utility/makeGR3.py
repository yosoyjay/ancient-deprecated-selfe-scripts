"""
Creates required set of *.gr3 files using standard values for nested
grids in the estuary.  Drag is constant!
"""
#-------------------------------------------------------------------------------
# Imports
#-------------------------------------------------------------------------------
from editDepthGR3 import change_gr3
import subprocess
import os

#-------------------------------------------------------------------------------
# Constants
#-------------------------------------------------------------------------------
BASEFILE = 'hgrid.gr3'
FILES = {'albedo.gr3'  : 0.06, 
         'diffmax.gr3' : 1.00, 
         'diffmin.gr3' : 1.0e-6,
         'drag.gr3'    : 0.002,
         'estuary.gr3' : 1.0,
				 'interpol.gr3': 2.0,
         'fluxflag.gr3': 1.0,
         's_nudge.gr3' : 0.0,
         't_nudge.gr3' : 0.0,
         'watertype.gr3' : 7.0,
         'xlsc.gr3'    : 0.5}

LL_BIN = '/home/workspace/users/lopezj/bin/spcs2ll'

WIND = 'windrot_geo2proj.gr3' 
WIND_BIN = '/home/workspace/users/lopezj/bin/rotate_wind_spcs2ll'

#-------------------------------------------------------------------------------
# Functions 
#-------------------------------------------------------------------------------
def makeGR3Files(dir, baseFile=BASEFILE):
  """Creates a set of *.gr3 files based on commonly used values in estuary only!

  Args:
    dir -- String of path to directory with the hgrid.gr3 and hgrid.ll files
    basefile -- String of path to hgrid.gr3 file that other files are based on
  Outputs:
    Required *.gr3 files for a model run.
  """
  old_wd = os.getcwd()
  os.chdir(dir)

  for f in FILES:
    change_gr3(baseFile, f, 100, FILES[f], FILES[f])
  
  args = [LL_BIN, '-input', 'hgrid.gr3', '-output', 'hgrid.ll', '-spcs2ll']
  proc = subprocess.Popen(args, shell=False)
  proc.wait()

  args = [WIND_BIN, '-input', 'hgrid.ll', '-output', 'windrot_geo2proj.gr3', '-ll2spcs']
  proc = subprocess.Popen(args, shell=False)
  proc.wait()

  os.chdir(old_wd)

#-------------------------------------------------------------------------------
# Main 
#-------------------------------------------------------------------------------
if __name__ == '__main__':
  import sys
  if len(sys.argv) == 1:
		raise Exception('Must provide path to run directory.')
  if len(sys.argv) == 2:
    makeGR3Files(sys.argv[1])
  if len(sys.argv) == 3:
    makeGR3Files(sys.argv[1], sys.argv[2])
