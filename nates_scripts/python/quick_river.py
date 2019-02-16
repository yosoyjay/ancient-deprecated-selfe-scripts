from numpy import *
from gr import *
from ob import *
import shape_procs
from model_util import *
import pdb
from model_reader import mreader
import meccs_file_util
import sys
from pylab import save, load
import meccs
import os

estuaryid  = "yaqalsea"
fin = "/home/hyde/models/yaqalsea/river_flow_4443_4495.dat"
out_base = "/home/hyde/models/yaqalsea/"

# set the time (which also dictates output file name)
sp = shape_procs.shape_procs()

rvr = sp.read_cs_data(fin)
mc = meccs.meccs()

#mean_river_flux
mean_flux = mc.get_time_avg_data(rvr['river_flow'][:,0], 900, 86400)
meccs_file_util.write_raw_data(out_base, "river_flux_"+str(4443)+"_"+str(53)+"_1.dat", \
                                 mean_flux, ["mean river flux"], "day:\t"+str(4443)+"\t"+str(53), 86400, "seconds", "m^3/s")                  
