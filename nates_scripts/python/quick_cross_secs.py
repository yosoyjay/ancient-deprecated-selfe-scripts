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

if len(sys.argv) != 9:
  sys.exit("usage: python quick_cross_secs.py [estuary id] [forecast directory] [shape_dir] [year] [day] [ndays] [data in dir] [out dir]")

estuaryid  = sys.argv[1]
base_dir   = sys.argv[2]
shape_dir = sys.argv[3]
yr         = int(sys.argv[4])
day        = int(sys.argv[5])
ndays      = int(sys.argv[6])
in_dir  = sys.argv[7]
out_base   = sys.argv[8]

# set the time (which also dictates output file name)
start_time = ydhms2corie(yr, day, 0, 15, 0)
end_time = start_time + ndays-0.25/24  #puts it right at the end of the day

sp = shape_procs.shape_procs()

# read all the shapes in: assume they are in the shape dir and named "points.txt", "cross_sections.txt", "transects.txt" and "regions.txt"

if os.access(shape_dir+"/cross_sections.txt", os.F_OK):  
  (cross_sections, cnames, crivers) = meccs_file_util.read_cs_file(shape_dir+"/cross_sections.txt", 0)
  do_cs=1
else:
  do_cs=0  

for i in arange(start_time, start_time+ndays, 1):
  if (i==start_time):
    rvr = sp.read_cs_data(in_dir+"/river_flow_"+str(int(i))+"_"+str(int(i+1))+".dat")
  else:
    tmp = sp.read_cs_data(in_dir+"/river_flow_"+str(int(i))+"_"+str(int(i+1))+".dat")
    rvr['river_flow'] = concatenate((rvr['river_flow'], tmp['river_flow']), 1)
    rvr['end_time'] = tmp['end_time']

if (do_cs):
  for i in arange(start_time, start_time+ndays, 1):
    if (i==start_time):
      csdata = sp.read_cs_data(in_dir+"/cross_section_data_"+str(int(i))+"_"+str(int(i+1))+".dat")
    else:
      tmp=sp.read_cs_data(in_dir+"/cross_section_data_"+str(int(i))+"_"+str(int(i+1))+".dat")
      for akey in tmp.keys():
        if (akey!='dt' and akey!='name' and akey!='end_time' and akey!='start_time'):        
          csdata[akey] = concatenate((csdata[akey], tmp[akey]), 1)
      csdata['end_time'] = tmp['end_time']

interval = 24.8 * 3600
tidal_period = 12.4 * 3600
dt = 900
mc = meccs.meccs()
ocean_salt = 33.5
ocean_den = 1024.0
salt_intr = 2.0 


#mean_river_flux
mean_flux = mc.get_time_avg_data(rvr['river_flow'][0,:], dt, interval)
meccs_file_util.write_raw_data(out_base, "mean_river_flux_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 mean_flux, ["mean river flux"], "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "m^3/s")                  

#cross_section calculations

cs_HnR = []
cs_froude = []
cs_flux = []
cs_area = []
cs_salt = []
cs_s_salt = []
cs_b_salt = []
cs_ac_surf_vel = []
cs_temp = []
cs_b_temp = []
cs_s_temp = []
cs_rho = []
cs_s_rho = []
cs_b_rho = []

for i in range(csdata['davg_salt'].shape[0]):
  pdb.set_trace()
  cs_HnR.append(mc.get_Hansen_Rattray(csdata['cs_surf_salt'][i,:], csdata['cs_bott_salt'][i,:], csdata['davg_salt'][i,:], \
                                      csdata['cs_surf_hvel_dir'][i,:], rvr['river_flow'][0,:]/csdata['cs_areas'][0,:], dt, interval))

  cs_froude.append(mc.get_time_avg_data(csdata['davg_hvel_mag'][i,:]/ \
                                        sqrt(9.8*(csdata['cs_bott_den'][i,:]-csdata['cs_surf_den'][i,:])/ocean_den*csdata['cs_depth'][i,:]), \
                                        dt, interval))
  cs_flux.append(csdata['hvel_flux'][i,:])
  cs_area.append(csdata['cs_areas'][i,:])
  cs_salt.append(csdata['davg_salt'][i,:])
  cs_s_salt.append(csdata['cs_surf_salt'][i,:])
  cs_b_salt.append(csdata['cs_bott_salt'][i,:])
  cs_ac_surf_vel.append(csdata['cs_surf_hvel_dir'][i,:])
  cs_temp.append(csdata['davg_temp'][i,:])
  cs_s_temp.append(csdata['cs_surf_temp'][i,:])
  cs_b_temp.append(csdata['cs_bott_temp'][i,:])
  cs_rho.append(csdata['davg_den'][i,:])
  cs_s_rho.append(csdata['cs_surf_den'][i,:])
  cs_b_rho.append(csdata['cs_bott_den'][i,:])  


tmp = zeros((cs_HnR[0].shape[1], len(cnames)*2), float)
tmpnames = []
for i in range(len(cnames)):
  tmp[:,i*2:i*2+2] = cs_HnR[i].T
  tmpnames.append("["+cnames[i]+": x, y]")
meccs_file_util.write_raw_data(out_base, "csection_Hansen_n_Rattray_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, tmpnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "dimensionless")

tmp = zeros((cs_froude[0].shape[0], len(cnames)), float)
for i in range(len(cnames)):
  tmp[:,i] = cs_froude[i]
meccs_file_util.write_raw_data(out_base, "csection_froude_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, cnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "dimensionless")

tmp = zeros((cs_flux[0].shape[0], len(cnames)), float)
for i in range(len(cnames)):
  tmp[:,i] = cs_flux[i]
meccs_file_util.write_raw_data(out_base, "csection_flux_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, cnames, "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "m^3/s")

tmp = zeros((cs_area[0].shape[0], len(cnames)), float)
for i in range(len(cnames)):
  tmp[:,i] = cs_area[i]
meccs_file_util.write_raw_data(out_base, "csection_area_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, cnames, "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "m^2")

tmp = zeros((cs_salt[0].shape[0], len(cnames)), float)
for i in range(len(cnames)):
  tmp[:,i] = cs_salt[i]
meccs_file_util.write_raw_data(out_base, "csection_avg_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, cnames, "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "psu")

tmp = zeros((cs_salt[0].shape[0], len(cnames)), float)
for i in range(len(cnames)):
  tmp[:,i] = cs_s_salt[i]
meccs_file_util.write_raw_data(out_base, "csection_surf_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, cnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((cs_salt[0].shape[0], len(cnames)), float)
for i in range(len(cnames)):
  tmp[:,i] = cs_b_salt[i]
meccs_file_util.write_raw_data(out_base, "csection_bott_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, cnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((cs_temp[0].shape[0], len(cnames)), float)
for i in range(len(cnames)):
  tmp[:,i] = cs_temp[i]
meccs_file_util.write_raw_data(out_base, "csection_depth_avg_temp_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, cnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((cs_temp[0].shape[0], len(cnames)), float)
for i in range(len(cnames)):
  tmp[:,i] = cs_s_temp[i]
meccs_file_util.write_raw_data(out_base, "csection_surf_temp_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, cnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((cs_temp[0].shape[0], len(cnames)), float)
for i in range(len(cnames)):
  tmp[:,i] = cs_b_temp[i]
meccs_file_util.write_raw_data(out_base, "csection_bott_temp_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, cnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((cs_salt[0].shape[0], len(cnames)), float)
for i in range(len(cnames)):
  tmp[:,i] = cs_rho[i]
meccs_file_util.write_raw_data(out_base, "csection_depth_avg_den_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, cnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((cs_salt[0].shape[0], len(cnames)), float)
for i in range(len(cnames)):
  tmp[:,i] = cs_s_rho[i]
meccs_file_util.write_raw_data(out_base, "csection_surf_den_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, cnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((cs_salt[0].shape[0], len(cnames)), float)
for i in range(len(cnames)):
  tmp[:,i] = cs_b_rho[i]
meccs_file_util.write_raw_data(out_base, "csection_bott_den_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, cnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

