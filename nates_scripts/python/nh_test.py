from numpy import *
from gr import *
from ob import *
from shape_procs import *
import pdb
from model_reader import mreader
import meccs_file_util

#agr = Gr()
#agr.readHGrid("/home/workspace/ccalmr/hindcasts/2004-01-14/run/hgrid.gr3")

pts = array([[334671.563016, 293717.116303], [338718.632401, 293613.689856], [342722.322867, 289452.065141], [343331.981229, 292526.175657], [349884.194234, 291364.688283], [350787.947752, 286265.475421], [358187.429680, 294594.189762], [358102.702788, 288361.818487]])

start_time = 2923.0104
end_time = 2923+(900./86400.)*96.*2



#ahdr = sz_header("/home/workspace/ccalmr37/hyde/coos/2004-01-11/run/1_elev.61")
#ahdr = sz_header("/home/workspace/ccalmr/hindcasts/2004-01-14/run/1_elev.61")
#ahdr.printHeader()
#ahdr.printVgrid()

#mr = mreader("/home/workspace/ccalmr/hindcasts/", bp=pts, start=2923.0104, end=2923+(900./86400.)*3., dvg=1, db=14)
#vel_dat = mr.get_model_data(dtype="hvel.64")
#salt_dat = mr.get_model_data(dtype="salt.63")

#(cs_dat, csnames) = meccs_file_util.read_cs_file("/home/users/hyde/columbia/col_meccs_cs_4.cs", 1)
#cs_dat = meccs_file_util.read_cs_file("/home/users/hyde/columbia/col_meccs_cs_4b.cs", 1)
#pts = meccs_file_util.read_pt_file("/home/users/hyde/columbia/col_meccs_4.pts", 0)

(cs_dat, csnames) = meccs_file_util.read_cs_file("/home/hyde/public_html/ourestuaries/shapes_data/hyde/yaqalsea/cross_sections.txt", 1)

meccs_file_util.write_cs_file("/home/hyde/public_html/ourestuaries/shapes_data/hyde/yaqalsea/test_cross_sections.txt", cs_dat, csnames)

(cs_dbg, dbg_names) = meccs_file_util.read_cs_file("/home/hyde/public_html/ourestuaries/shapes_data/hyde/yaqalsea/test_cross_sections.txt", 1)

#sp = shape_procs()
#for i in range(len(areas)):
#  sp.add_shape(areas[i], "volume", 0, "vol_"+str(i))
#data = sp.gen_volume_data("/home/workspace/ccalmr/hindcasts/", 2923.0104, 2923+(900./86400.)*2., 16, "base", 1)

sp = shape_procs()

for i in range(len(cs_dat)):
  sp.add_shape(cs_dat[i], "cross_sec", 50, "cs"+str(i))
  
data = sp.gen_cross_section_data("/home/workspace/ccalmr/hindcasts/", start_time, end_time, 16, "base", 1)

data["start_time"] = start_time
data["end_time"] = end_time
sp.save_cs_data(data, "test_save2.sp")

#data2 = sp.read_cs_data("test_save.sp")


#for i in range(pts.shape[0]):
#  sp.add_shape(pts[i], "point", 0, "pt"+str(i))

#data = sp.gen_point_data("/home/workspace/ccalmr/hindcasts/", 2923.0104, 2923+(900./86400.)*2., 16, "base", 1)

dbg = 1
