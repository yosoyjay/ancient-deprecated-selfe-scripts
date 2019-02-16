import sys
sys.path.append('/home/hyde/ourestuaries/pythonlib/cmop')
from numpy import *
from gr import *
from ob import *
import shape_procs
from model_util import *
#import pdb
from model_reader import mreader
import meccs_file_util
#import sys
from pylab import save, load
import os.path
import pdb

if len(sys.argv) < 6:
  sys.exit("usage: python quick_rivers.py [estuary id] [processed data directory] [year] [day] [ndays] [out dir]")

estuaryid  = sys.argv[1]
base_dir   = sys.argv[2]
yr         = int(sys.argv[3])
day        = int(sys.argv[4])
ndays      = int(sys.argv[5])
out_base   = sys.argv[6]

# set the time (which also dictates output file name)
start_time = ydhms2corie(yr, day, 0, 15, 0)
end_time = start_time + ndays-0.25/24  #puts it right at the end of the day

print "start_time: "+str(start_time)

# get the river flow data
sp = shape_procs.shape_procs()
mr = mreader(base_dir+"/"+estuaryid+"/", array([]), start_time, end_time, 0, 16, 1, 0)
river_flow=mr.get_river_flux(base_dir+"/"+estuaryid+"/", start_time, end_time, 16, 900, 0)
data = {}
data["river_flow"] = river_flow
data["start_time"] = start_time
data["end_time"] = end_time
data["dt"] = 900
sp.save_cs_data(data, out_base+"/river_flow_"+str(int(start_time))+"_"+str(int(end_time))+".dat")
