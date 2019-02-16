import sys
sys.path.append('/home/hyde/ourestuaries/pythonlib/cmop')
from numpy import *
from gr import *
from ob import *
import shape_procs
from model_util import *
from model_reader import mreader
import meccs_file_util
import meccs

if len(sys.argv) < 5:
  sys.exit("usage: python get_net_flux.py [processed flux file] [out dir] [name base] [time avg period (seconds)] [river index]...")

flux_file = sys.argv[1]
out_dir = sys.argv[2]
name_base = sys.argv[3]
interval = int(sys.argv[4])

rivers = []
for i in range(5,len(sys.argv)):
  rivers.append(int(sys.argv[i]))

fin = flux_file

# set the time (which also dictates output file name)
sp = shape_procs.shape_procs()

rvr = sp.read_cs_data(fin)
mc = meccs.meccs()

#mean_river_flux
for i in rivers:
  if i==rivers[0]:
    mean_flux = mc.get_time_avg_data(rvr['river_flow'][:,i], 900, interval)
  else:
    mean_flux = mean_flux+mc.get_time_avg_data(rvr['river_flow'][:,i], 900, interval)

meccs_file_util.write_raw_data(out_dir, name_base+"_river_flux.dat", \
                                 mean_flux, ["mean river flux"], " ", interval, "seconds", "m^3/s")                  
