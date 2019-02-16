import sys
sys.path.append('/home/users/hyde/ourestuaries/pythonlib/cmop')
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
import interpolate

if len(sys.argv) != 8:
  sys.exit("usage: python do_daily_process.py [estuary id] [base_dir] [shape_dir] [year] [day] [ndays] [out dir]")

estuaryid  = sys.argv[1]
base_dir   = sys.argv[2]
shape_dir = sys.argv[3]
yr         = int(sys.argv[4])
day        = int(sys.argv[5])
ndays      = int(sys.argv[6])
out_dir   = sys.argv[7]

data_types = ["salt", "temp", "hvel_u", "hvel_v", "elev", "zcor"]

#for computing significant equality
small = 0.0001

# set the time (which also dictates output file name)
start_time = ydhms2corie(yr, day, 0, 15, 0)
end_time = start_time + ndays-0.25/24  #puts it right at the end of the day

#time intervals and other constants
interval = 24.8 * 3600
tidal_period = 12.4 * 3600
dt = 900
mc = meccs.meccs()
ocean_salt = 33.5
ocean_den = 1024.0
salt_intr = 2.0 
tran_space_interval = 10

#variables
salt_intrusion = []
tr_froude = []
surf_vel = []
da_vel = []
da_mag = []
surf_salt = []
bott_salt = []
surf_den = []
bott_den = []
davg_den = []
davg_salt = []
out_xy = []

sp = shape_procs.shape_procs()

(transects, tnames, trivers) = meccs_file_util.read_cs_file(shape_dir+"/transects.txt", 1)

tdata = []
for i in range(ndays):
  curr_time = start_time+i
  in_dir = "%s/%s/%d-%03d/process/" % (base_dir, estuaryid, yr, day+i)    
  for j in range(len(transects)):

    if (i==0):
      tdata.append(sp.read_cs_data(in_dir+"/"+tnames[j]+"_transect_data.dat"))
    else:
      tmp = sp.read_cs_data(in_dir+"/"+tnames[j]+"_transect_data.dat")
      for akey in tmp.keys():
        if (akey!='dt' and akey!='name' and akey!='end_time' and akey!='start_time'):
          tdata[j][akey] = concatenate((tdata[j][akey], tmp[akey]), 1)
      tdata[j]['end_time'] = tmp['end_time']


for i in range(len(transects)):
  sp.add_shape(transects[i], "transect", 50, tnames[i])
  temp_shape = sp.get_shape(tnames[i])

  salt_intrusion.append(mc.get_salinity_intrusion(temp_shape, tdata[i]['dmax_salt'], dt, interval, salt_intr))
  out_xy.append(temp_shape[0::tran_space_interval, 0:2])
  surf_vel.append(tdata[i]['surf_hvel_dir'][0::tran_space_interval,:])
  da_vel.append(tdata[i]['davg_hvel_dir'][0::tran_space_interval,:])
  da_mag.append(tdata[i]['davg_hvel_mag'][0::tran_space_interval,:])
  surf_salt.append(tdata[i]['surf_salt'][0::tran_space_interval,:])
  bott_salt.append(tdata[i]['bott_salt'][0::tran_space_interval,:])
  davg_salt.append(tdata[i]['davg_salt'][0::tran_space_interval,:])
  surf_den.append(tdata[i]['surf_den'][0::tran_space_interval,:])
  bott_den.append(tdata[i]['bott_den'][0::tran_space_interval,:])  
  davg_den.append(tdata[i]['davg_den'][0::tran_space_interval,:])

  tstep_int = floor(interval/dt)  #keeping it simple and slightly innacurate
  tmp = zeros((tdata[i]['bott_den'][0::tran_space_interval,:].shape[0], int(tdata[i]['bott_den'].shape[1]/tstep_int)), float)
  dpth = tdata[i]['depth'][0::tran_space_interval,:]

  for j in range(tdata[i]['bott_den'][0::tran_space_interval,:].shape[0]):
    tmp[j,:] = mc.get_time_avg_data(da_mag[i][j,:]/sqrt(9.8*(bott_den[i][j,:]-surf_den[i][j,:])/1024*dpth[j,:]), dt, interval)
  tr_froude.append(tmp)

# output transect data:

for i in range(len(transects)):
  meccs_file_util.write_raw_data(out_dir, "transect_"+tnames[i]+"_salt_intrusion_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 salt_intrusion[i].T, ["salinity intrusion"], "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "meters")
  meccs_file_util.write_raw_data(out_dir, "transect_"+tnames[i]+"_output_xy_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 out_xy[i], ["x", "y"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "N/A", "meters")
  meccs_file_util.write_raw_data(out_dir, "transect_"+tnames[i]+"_surf_vel_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 surf_vel[i].T, ["surface velocity"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "m/s")
  meccs_file_util.write_raw_data(out_dir, "transect_"+tnames[i]+"_ac_davg_vel_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 da_vel[i].T, ["along channel depth averaged velocity"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "m/s")
  meccs_file_util.write_raw_data(out_dir, "transect_"+tnames[i]+"_mag_depth_avg_vel_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 da_mag[i].T, ["magnitude of depth averaged velocity"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "m/s")
  meccs_file_util.write_raw_data(out_dir, "transect_"+tnames[i]+"_surface_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 surf_salt[i].T, ["surface salinity"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "psu")
  meccs_file_util.write_raw_data(out_dir, "transect_"+tnames[i]+"_bottom_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 bott_salt[i].T, ["bottom salinity"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "psu")
  meccs_file_util.write_raw_data(out_dir, "transect_"+tnames[i]+"_davg_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 davg_salt[i].T, ["depth averaged salinity"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "psu")
  meccs_file_util.write_raw_data(out_dir, "transect_"+tnames[i]+"_surface_rho_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 surf_den[i].T, ["surface density"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "kg/m^3")
  meccs_file_util.write_raw_data(out_dir, "transect_"+tnames[i]+"_bottom_rho_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 bott_den[i].T, ["bottom density"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "kg/m^3")
  meccs_file_util.write_raw_data(out_dir, "transect_"+tnames[i]+"_bottom_rho_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 bott_den[i].T, ["bottom density"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "kg/m^3")
  meccs_file_util.write_raw_data(out_dir, "transect_"+tnames[i]+"_davg_rho_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 davg_den[i].T, ["depth averaged density"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "kg/m^3")
  meccs_file_util.write_raw_data(out_dir, "transect_"+tnames[i]+"_froude_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tr_froude[i].T, ["interfacial Froude #"], "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "kg/m^3")

