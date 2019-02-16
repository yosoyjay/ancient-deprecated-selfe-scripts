from numpy import *
from gr import *
from ob import *
from shape_procs import *
from model_util import *
import pdb
from model_reader import mreader
import meccs_file_util
import sys
import re

# command line arguments, in order:
#  - esutuary id
#  - model base directory
#  - model directory format
#  - start year
#  - start day
#  - end year
#  - end day
#  - shape type
#  - shape file

if len(sys.argv) != 10:
  sys.exit("usage: python do_get_shapes.py [estuary id] [model base directory] [model directory format - e.g. yyyy-ww-17 for hindcast, db 17] [start year] [start day] [end year] [end day] [shape type] [shape file]")

estuaryid  = sys.argv[1]
base_dir   = sys.argv[2]
dir_format = sys.argv[3]
syr        = sys.argv[4]
sday       = sys.argv[5]
eyr        = sys.argv[6]
eday       = sys.argv[7]
shape_type = sys.argv[8]
shape_file = sys.argv[9]

start_time = ydhms2corie(int(syr), int(sday), 0, 15, 0)
end_time = ydhms2corie(int(eyr), int(eday), 0, 0, 0)

hindcast = 1
ndaysDir = 7
if re.match('yyyy-ddd', dir_format):
  hindcast = 0
  ndaysDir = 1

#print "start_time: "+str(start_time)+"\n"

if shape_type == "point":
  pass
elif shape_type == "cross_section":
  cs_dat = meccs_file_util.read_cs_file(shape_file, 1)

  do_cs_norms = 0
  sp = shape_procs()
  for i in range(len(cs_dat)):
    sp.add_shape(cs_dat[i], "cross_sec", 0, "cs"+str(i))
    tmp = cs_dat[i]
    if tmp[0][2]==0 and tmp[0][3]==0:    #assume if the first normal vectors are both 0, then the norms need to be calculated
      do_cs_norms=1

  if do_cs_norms:
    sp.gen_channel_norms(base_dir, start_time, start_time+2, "cross_sec", 16, ndaysDir, hindcast)

  data = sp.gen_cross_section_data(base_dir, start_time, end_time, 16, "base", 1, ndaysDir, hindcast)

  data["start_time"] = start_time
  data["end_time"] = end_time
  sp.save_cs_data(data, estuaryid+"_"+str(syr)+"_"+str(sday)+"_to_"+str(eyr)+"_"+str(eday)+".sp");

#  data = sp.read_cs_data(estuarid+estuaryid+"_"+str(syr)+"_"+str(sday)+"_to_"+str(eyr)+"_"+str(eday)+".sp");
  
elif shape_type == "transect":
  ts_dat = meccs_file_util.read_cs_file(shape_file, 1)  #transect and cross section have same file format
  
  do_norms = 0
  sp = shape_procs()
  for i in range(len(ts_dat)):
    sp.add_shape(ts_dat[i], shape_type, 0, "ts"+str(i))
    tmp = ts_dat[i]
    if tmp[0][2]==0 and tmp[0][3]==0:    #assume if the first normal vectors are both 0, then the norms need to be calculated
      do_norms = 1
  if do_norms:
    sp.gen_channel_norms(base_dir, start_time, start_time+2, "transect", 16, ndaysDir, hindcast)

  for i in range(len(ts_dat)):  
    data = sp.gen_transect_data(base_dir, start_time, end_time, "ts"+str(i), 16, 1, ndaysDir, hindcast)

    data["start_time"] = start_time
    data["end_time"] = end_time
    sp.save_cs_data(data, estuaryid+"_ts"+str(i)+"_"+str(syr)+"_"+str(sday)+"_to_"+str(eyr)+"_"+str(eday)+".sp");
  
elif shape_type == "volume":

  area_dat = meccs_file_util.read_area_file(shape_file)

  sp = shape_procs()
  for i in range(len(area_dat)):
    sp.add_shape(area_dat[i], "volume", 0, "vol"+str(i))

  data = sp.gen_volume_data(base_dir, start_time, end_time, 16, "base", 1, 1)

  pdb.set_trace()
  
  data["start_time"] = start_time
  data["end_time"] = end_time
  sp.save_cs_data(data, estuaryid+"_volume_"+str(syr)+"_"+str(sday)+"_to_"+str(eyr)+"_"+str(eday)+".sp");

else:
  sys.exit("error, invalid shape type: "+shape_type)



