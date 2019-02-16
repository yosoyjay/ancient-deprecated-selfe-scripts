from numpy import *
from gr import *
from ob import *
import shape_procs
from model_util import *
import pdb
from model_reader import mreader
import meccs_file_util
import sys
import meccs
from pylab import save, load

if len(sys.argv) != 10:
  sys.exit("usage: python do_class_num_test.py [estuary id] [processed data directory] [start year] [start day] [end year] [end day] [shape type] [shape file]")

estuaryid  = sys.argv[1]
base_dir   = sys.argv[2]
syr        = sys.argv[3]
sday       = sys.argv[4]
eyr        = sys.argv[5]
eday       = sys.argv[6]
ctype      = sys.argv[7]
sfile      = sys.argv[8]
out_base   = sys.argv[9]

if ctype == "intrusion":
  s_dat = meccs_file_util.read_cs_file(sfile, 1)

  mec = meccs.meccs()
  sp = shape_procs.shape_procs()
  
  intrusion = []
  for i in range(len(s_dat)):
    data = sp.read_cs_data(base_dir+"/"+estuaryid+"_ts"+str(i)+"_"+str(syr)+"_"+str(sday)+"_to_"+str(eyr)+"_"+str(eday)+".sp")
    intrusion.append(mec.get_salinity_intrusion(s_dat[i], data["dmax_salt"], data["dt"], 24.8 * 3600, 2))

  pdb.set_trace()
  dbg=1

if ctype == "hansen_rattray":
  s_dat = meccs_file_util.read_cs_file(sfile, 1)

  mec = meccs.meccs()
  sp = shape_procs.shape_procs()

  data = sp.read_cs_data(base_dir+"/"+estuaryid+"_"+str(syr)+"_"+str(sday)+"_to_"+str(eyr)+"_"+str(eday)+".sp")

  su_u = zeros(data["cs_surf_salt"].shape, float)
  ds_s = zeros(data["cs_surf_salt"].shape, float)  
  legendstr = []

  for i in range(len(s_dat)):
    dat = mec.get_Hansen_Rattray(data["cs_surf_salt"][i,:], data["cs_bott_salt"][i,:], data["davg_salt"][i,:], \
                             data["cs_surf_hvel_dir"][i,:], data["davg_hvel_mag"][i,:], data["dt"], 24.8 * 3600)

    save(out_base+data["names"][i][0]+".hr", dat, fmt="%20.10f")
    
    legendstr.append(data["names"][i][0])
    su_u[i,:] = dat[0,:]
    ds_s[i,:] = dat[1,:]    

  mec.hnr_plot(su_u, ds_s, "test HNR plot", legendstr)  

if ctype == "river_flux":
  start_time = ydhms2corie(int(syr), int(sday), 0, 15, 0)
  end_time = ydhms2corie(int(eyr), int(eday), 0, 0, 0)

  mr = mreader(base_dir, array([]), start_time, end_time)
  data = mr.get_river_flux(base_dir, start_time, end_time, 11, 900, [0, 1], 1)  
