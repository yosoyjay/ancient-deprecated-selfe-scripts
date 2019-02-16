from numpy import *
from gr import *
from ob import *
import shape_procs
from model_util import *
import pdb
from model_reader import mreader
import meccs_file_util
import sys
import meccs
from pylab import save, load

if len(sys.argv) != 9:
  sys.exit("usage: python do_class_num_test.py [estuary id] [processed data directory] [shape_dir] [year] [day] [ndays] [shape file] [out dir]")

estuaryid  = sys.argv[1]
base_dir   = sys.argv[2]
shape_dir  = sys.argv[3]
yr         = sys.argv[4]
day        = sys.argv[5]
ndays      = sys.argv[6]
sfile      = sys.argv[7]
out_base   = sys.argv[8]

# read all the shapes in: assume they are in the shape dir and named "points.txt", "cross_sections.txt", "transects.txt" and "regions.txt"
points = read_pt_file

# set the time (which also dictates output file name)
start_time = ydhms2corie(int(syr), int(sday), 0, 15, 0)
end_time = ydhms2corie(int(eyr), int(eday), 0, 0, 0)

# Flow Ratio


# Tidal Exchange


# Filling Time


# Fresh Water Flushing Time


# Fraction of Saltwater


# Intrusion Length


# Interfacial Froude #


# Cross section flux vs. elev area


# Interfacial Froude (cs & point)


# Prandle variables (cs & point)


# Hansen & Rattray (cs & point)


# volume, salt volume and heat volume time-series
