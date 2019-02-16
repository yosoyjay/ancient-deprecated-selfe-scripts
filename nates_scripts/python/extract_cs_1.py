import sys
sys.path.append('/home/users/hyde/ourestuaries/pythonlib/cmop')
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

if len(sys.argv) < 9:
  sys.exit("usage: python extract_cs_1.py [estuary id] [processed data directory] [shape_dir] [year] [day] [ndays] [out dir] [full refresh]")

estuaryid  = sys.argv[1]
base_dir   = sys.argv[2]
shape_dir  = sys.argv[3]
yr         = int(sys.argv[4])
day        = int(sys.argv[5])
ndays      = int(sys.argv[6])
out_base   = sys.argv[7]
do_refresh = sys.argv[8]

# read all the shapes in: assume they are in the shape dir and named "points.txt", "cross_sections.txt", "transects.txt" and "regions.txt"
#(points, pnames) = meccs_file_util.read_cs_file(shape_dir+"/points.txt", 0)   #uses cs file format
#(cross_sections, cnames) = meccs_file_util.read_cs_file(shape_dir+"/cross_sections.txt", 0)
#(regions, rnames) = meccs_file_util.read_cs_file(shape_dir+"/regions.txt", 0)
#(transects, tnames) = meccs_file_util.read_cs_file(shape_dir+"/transects.txt", 0)

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

# Cross sections data: flux, area, cross-sectionally averaged and normalized surface and bottom for temp, salt density and hvel (surf onl
if (os.path.exists(shape_dir+"/cross_sections.txt")):
  sp = shape_procs.shape_procs()

  (cross_sections, cnames, crivers) = meccs_file_util.read_cs_file(shape_dir+"/cross_sections.txt", 1)
  if len(cross_sections)>0:
    for i in range(len(cross_sections)):
      sp.add_shape(cross_sections[i], "cross_sec", 50, cnames[i])

    if (meccs_file_util.file_exists(shape_dir+"/cross_sections.ob") and int(do_refresh)==0):
      print "using cross section ob file "+shape_dir+"cross_sections.ob\n"      
      sp.ob_file=shape_dir+"/cross_sections.ob"
    else:
      print "calculating new cross section ob file\n"
      agr = Gr()
      [year, day, ts] = model_util.corie2ydn(start_time)
      dname = "%s/%d-%03d/run/" % (base_dir+"/"+estuaryid+"/", year, day)      
      agr.readHGrid(dname+"hgrid.gr3")
      aob = Ob()
      sppts = sp.get_pts("cross_sec")
      
      aob.grid_init(agr, array([sppts[0], sppts[1]]).T)
      aob.save_ob(shape_dir+"/cross_sections.ob")
      sp.ob_file = shape_dir+"/cross_sections.ob"

    if (mean(abs(cross_sections[0][:,2]))==0 and mean(abs(cross_sections[0][:,3]))==0) or int(do_refresh)==1:
      sp.gen_channel_norms(base_dir+"/"+estuaryid+"/", start_time, start_time+2, "cross_sec", 16, 1, 0)
      meccs_file_util.write_cs_file(shape_dir+"/cross_sections.txt", sp.shapes, cnames, crivers)

    data = sp.gen_cross_section_data(base_dir+"/"+estuaryid+"/", start_time, end_time, 16, "base", 1, 1, 0)
    
    data["start_time"] = start_time
    data["end_time"] = end_time
    sp.save_cs_data(data, out_base+"/cross_section_data_"+str(int(start_time))+"_"+str(int(end_time))+".dat");
    
