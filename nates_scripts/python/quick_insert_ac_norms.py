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
import ourestuaries_db

shape_dir = '/home/workspace/ccalmr38/hindcasts/yaqalsea_files/'
model_dir = '/home/workspace/local0/forecasts/yaqalsea/'
year = 2008
sday = 95
ndays = 5
start_time = ydhms2corie(year, sday, 0, 15, 0)
end_time = start_time + ndays-0.25/24  #puts it right at the end of the day

(transects, tnames, trivers) = meccs_file_util.read_cs_file(shape_dir+"/transects.txt", 1)
(points, pnames, privers) = meccs_file_util.read_cs_file(shape_dir+"/points.txt", 1)
(csecs, cnames, crivers) = meccs_file_util.read_cs_file(shape_dir+"/cross_sections.txt", 1)

oe = ourestuaries_db.oe_db()
oe.connect("cdb01.stccmop.org", "ourestuaries", "test_hyde", "poop77")

#(names, pts, norms, rivers) = oe.get_points_by_eid('yaqalsea')
#(names, runids, unames) = oe.get_shape_info_by_eid('yaqalsea', 'cross_section')
#points = 
#for i in range(len(names)):
  
  
#pdb.set_trace()

#sp = shape_procs.shape_procs()
#for i in range(pts.shape[0]):
#  sp.add_shape(array([hstack((pts[i],[0,0]))])  , "point", 0, names[i])

pdb.set_trace()
  
sp.gen_channel_norms(model_dir, start_time, end_time, "point", 16, 1, 0)

meccs_file_util.write_cs_file(shape_dir+"/temp_pts.txt", sp.shapes, names)

#(pts, norms, rivers, eid) = oe.get_trans_or_cs("ts1", "yaqadj", "hyde", "transect")
#oe.set_transect(transects[0][:,0:2], transects[0][:,2:4], "tstest", rivers, eid, "yaqadj", "hyde")

dbg=1


