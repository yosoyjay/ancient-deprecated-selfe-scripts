import sys
sys.path.append('/home/hyde/ourestuaries/pythonlib/cmop')
from numpy import *
from gr import *
from ob import *
import shape_procs
from model_util import *
from model_reader import mreader
import meccs_file_util
from pylab import save, load
import os.path
import pdb
import 

if len(sys.argv) < 10:
  sys.exit("usage: python setup_hindcast_extraction [base_dir] [run name] [target_dir] [do points 0 or 1] [do transects 0 or 1] [do cross sections 0 or 1] [do regions 0 or 1] [year - for calculating along-channel normal vectors] [start day (normal vectors)] [end day (normal vectors)]")

base_dir = sys.argv[1]
run_name = sys.argv[2]
target_dir = sys.argv[3]
do_points = integer(sys.argv[4])
do_trans = integer(sys.argv[5])
do_csecs = integer(sys.argv[6])
do_regs = integer(sys.argv[7])
year = integer(sys.argv[8])
sday = integer(sys.argv[9])
eday = integer(sys.argv[10])

# todo:
#   - for each do_[shape] == 1
#      - query database for shape type and run_name
#      - generate along-channel normal vectors (except for regions)
#      - write the result out using meccs_file_util.write_cs_file
#      - later, add upstream_rivers array to the shapes and include it in the write_cs_file call


