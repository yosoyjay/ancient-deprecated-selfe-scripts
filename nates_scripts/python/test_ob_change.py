import sys
sys.path.append('/home/users/hyde/ourestuaries/pythonlib/cmop')
from numpy import *
from gr import *
from ob import *
import shape_procs
import shape_procs_test
from model_util import *
#import pdb
from model_reader import mreader
import meccs_file_util
#import sys
from pylab import save, load
import os.path
import pdb
import time

if len(sys.argv) < 8:
  sys.exit("usage: python test_ob_change.py [estuary id] [processed data directory] [shape_dir] [year] [day] [ndays] [out dir]")

estuaryid  = sys.argv[1]
base_dir   = sys.argv[2]
shape_dir  = sys.argv[3]
yr         = int(sys.argv[4])
day        = int(sys.argv[5])
ndays      = int(sys.argv[6])
out_base   = sys.argv[7]

line_res = 100;

# set the time (which also dictates output file name)
start_time = ydhms2corie(yr, day, 0, 15, 0)
end_time = start_time + ndays-0.25/24  #puts it right at the end of the day

print "start_time: "+str(start_time)

# Cross sections data: flux, area, cross-sectionally averaged and normalized surface and bottom for temp, salt density and hvel (surf onl

(cross_sections, cnames, crivers) = meccs_file_util.read_cs_file(shape_dir+"/cross_sections.txt", 1)
 
sp = shape_procs.shape_procs()

for i in range(len(cross_sections)):
  sp.add_shape(cross_sections[i], "cross_sec", line_res, cnames[i])

t1 = time.time()

data = sp.gen_cross_section_data(base_dir+"/"+estuaryid+"/", start_time, end_time, 16, "base", 1, 1, 0)

t2 = time.time()

print "time for old style: "+str(t2-t1)+"\n";

data["start_time"] = start_time
data["end_time"] = end_time
sp.save_cs_data(data, out_base+"/test_data_1.dat");
    
# using the new calculation style

sp = shape_procs_test.shape_procs_test()

for i in range(len(cross_sections)):
  sp.add_shape(cross_sections[i], "cross_sec", line_res, cnames[i])

agr = Gr()
[year, day, ts] = model_util.corie2ydn(start_time)
dname = "%s/%d-%03d/run/" % (base_dir+"/"+estuaryid+"/", year, day)      
agr.readHGrid(dname+"hgrid.gr3")
aob = Ob()
sppts = sp.get_pts("cross_sec")

aob.grid_init(agr, array([sppts[0], sppts[1]]).T)
aob.save_ob(shape_dir+"/cross_sections.ob")
sp.ob_file = shape_dir+"/cross_sections.ob"

t1 = time.time()

data = sp.gen_cross_section_data(base_dir+"/"+estuaryid+"/", start_time, end_time, 16, "base", 1, 1, 0)

t2 = time.time()

print "time for new style: "+str(t2-t1)+"\n";

data["start_time"] = start_time
data["end_time"] = end_time
sp.save_cs_data(data, out_base+"/test_data_3.dat");

# test loading and unloading
#
#sp = shape_procs_test.shape_procs_test()
#
#for i in range(len(cross_sections)):
#  sp.add_shape(cross_sections[i], "cross_sec", line_res, cnames[i])
#
#t1 = time.time()
#
#data = sp.gen_cross_section_data(base_dir+"/"+estuaryid+"/", start_time, end_time, 16, "base", 1, 1, 0)
#
#t2 = time.time()
#
#print "time for new style: "+str(t2-t1)+"\n";
#
#data["start_time"] = start_time
#data["end_time"] = end_time
#sp.save_cs_data(data, out_base+"/test_data_2.dat");
