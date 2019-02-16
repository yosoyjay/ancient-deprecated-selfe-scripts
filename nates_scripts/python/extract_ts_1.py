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

if len(sys.argv) < 8:
  sys.exit("usage: python extract_ts_1.py [forecast directory] [shape_dir] [year] [day] [ndays] [out dir] [full refresh]")

base_dir   = sys.argv[1]
shape_dir  = sys.argv[2]
yr         = int(sys.argv[3])
day        = int(sys.argv[4])
ndays      = int(sys.argv[5])
out_base   = sys.argv[6]
do_refresh = sys.argv[7]

# set the time (which also dictates output file name)
start_time = ydhms2corie(yr, day, 0, 15, 0)
end_time = start_time + ndays-0.25/24  #puts it right at the end of the day

# get the river flow data
sp = shape_procs.shape_procs()
mr = mreader(base_dir+"/", array([]), start_time, end_time, 0, 16, 1, 0)
river_flow=mr.get_river_flux(base_dir+"/", start_time, end_time, 16, 900, 0)
data = {}
data["river_flow"] = river_flow
data["start_time"] = start_time
data["end_time"] = end_time
data["dt"] = 900
sp.save_cs_data(data, out_base+"/river_flow_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat")

# transect calculations: min, max salinity, surface salt, bottom salt, depth-avg salt,
# along-channel min, max, depth_avg velocity by dir & magnitude

if (os.path.exists(shape_dir+"/transects.txt")):
  sp = shape_procs.shape_procs()
  
  (transects, tnames, trivers) = meccs_file_util.read_cs_file(shape_dir+"/transects.txt", 0)
  if len(transects)>0:
    for i in range(len(transects)):
      sp.add_shape(transects[i], "transect", 50, tnames[i])

    if (meccs_file_util.file_exists(shape_dir+"/transects.ob") and int(do_refresh)==0):
#      print "using transect ob file "+shape_dir+"transects.ob\n"      
      sp.ob_file=shape_dir+"/transects.ob"
    else:
#      print "creating new transect ob file\n"
      agr = Gr()
      [year, day, ts] = model_util.corie2ydn(start_time)
      dname = "%s/%d-%03d/run/" % (base_dir+"/", year, day)      
      agr.readHGrid(dname+"hgrid.gr3")
      aob = Ob()
      sppts = sp.get_pts("transect")
      aob.grid_init(agr, array([sppts[0], sppts[1]]).T)
      aob.save_ob(shape_dir+"/transects.ob")
      sp.ob_file = shape_dir+"/transects.ob"

    if (mean(abs(transects[0][:,2]))==0 and mean(abs(transects[0][:,3]))==0) or int(do_refresh)==1:
      sp.gen_channel_norms(base_dir+"/", start_time, start_time+2, "transect", 16, 1, 0)
      meccs_file_util.write_cs_file(shape_dir+"/transects.txt", sp.shapes, tnames, trivers)      
    
    for i in range(len(transects)):  
      data = sp.gen_transect_data(base_dir+"/", start_time, end_time, tnames[i], 16, 1, 1, 0)
    
      data["start_time"] = start_time
      data["end_time"] = end_time
      sp.save_cs_data(data, out_base+"/transect_"+tnames[i]+"_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat");
