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
  sys.exit("usage: python cnc_set_1.py [estuary id] [forecast directory] [shape_dir] [year] [day] [ndays] [data in dir] [out dir]")

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
if os.access(shape_dir+"/points.txt", os.F_OK):
  (points, pnames, privers) = meccs_file_util.read_cs_file(shape_dir+"/points.txt", 0)   #uses cs file format
  do_points=1
else:
  do_points=0
if os.access(shape_dir+"/cross_sections.txt", os.F_OK):  
  (cross_sections, cnames, crivers) = meccs_file_util.read_cs_file(shape_dir+"/cross_sections.txt", 0)
  do_cs=1
else:
  do_cs=0  
if os.access(shape_dir+"/transects.txt", os.F_OK):  
  (transects, tnames, trivers) = meccs_file_util.read_cs_file(shape_dir+"/transects.txt", 0)
  do_transects=1
else:
  do_transects=0  
if os.access(shape_dir+"/regions.txt", os.F_OK):  
  (regions, rnames, rrivers) = meccs_file_util.read_cs_file(shape_dir+"/regions.txt", 0)
  do_regions=1
else:
  do_regions=0

for i in arange(start_time, start_time+ndays, 1):
  if (i==start_time):
    rvr = sp.read_cs_data(in_dir+"/river_flow_"+str(int(i))+"_"+str(int(i+1))+".dat")
  else:
    tmp = sp.read_cs_data(in_dir+"/river_flow_"+str(int(i))+"_"+str(int(i+1))+".dat")
    rvr['river_flow'] = concatenate((rvr['river_flow'], tmp['river_flow']), 1)
    rvr['end_time'] = tmp['end_time']

if (do_transects):
  tdata = []
  for i in arange(start_time, start_time+ndays, 1):
    for j in range(len(transects)):
      if (i==start_time):
        tdata.append(sp.read_cs_data(in_dir+"/transect_"+tnames[j]+"_"+str(int(i))+"_"+str(int(i+1))+".dat"))
      else:
        tmp = sp.read_cs_data(in_dir+"/transect_"+tnames[j]+"_"+str(int(i))+"_"+str(int(i+1))+".dat")
        for akey in tmp.keys():
          if (akey!='dt' and akey!='name' and akey!='end_time' and akey!='start_time'):
            tdata[j][akey] = concatenate((tdata[j][akey], tmp[akey]), 1)
        tdata[j]['end_time'] = tmp['end_time']

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

if (do_points):
  for i in arange(start_time, start_time+ndays, 1):
    if (i==start_time):  
      ptdata = sp.read_cs_data(in_dir+"/point_data_"+str(int(i))+"_"+str(int(i+1))+".dat")
    else:
      tmp = sp.read_cs_data(in_dir+"/point_data_"+str(int(i))+"_"+str(int(i+1))+".dat")
      for akey in tmp.keys():
        if (akey!='dt' and akey!='name' and akey!='end_time' and akey!='start_time'):        
          ptdata[akey] = concatenate((ptdata[akey], tmp[akey]), 1)
      ptdata['end_time'] = tmp['end_time']
      

if (do_regions):
  for i in arange(start_time, start_time+ndays, 1):
    if (i==start_time):  
      vdata = sp.read_cs_data(in_dir+"/volume_data_"+str(int(i))+"_"+str(int(i+1))+".dat")
    else:
      tmp = sp.read_cs_data(in_dir+"/volume_data_"+str(int(i))+"_"+str(int(i+1))+".dat")
      for akey in tmp.keys():
        if (akey!='dt' and akey!='name' and akey!='end_time' and akey!='start_time'):        
          vdata[akey] = concatenate((vdata[akey], tmp[akey]), 1)
      vdata['end_time'] = tmp['end_time']

interval = 24.8 * 3600
tidal_period = 12.4 * 3600
dt = 900
mc = meccs.meccs()
ocean_salt = 33.5
ocean_den = 1024.0
salt_intr = 2.0 


#transect data - do every 10 pts
tran_space_interval = 10
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

  tstep_int = ceil(interval/dt)  #keeping it simple and slightly innacurate
  tmp = zeros((tdata[i]['bott_den'][0::tran_space_interval,:].shape[0], int(tdata[i]['bott_den'].shape[1]/tstep_int)), float)
  dpth = tdata[i]['depth'][0::tran_space_interval,:]
  for j in range(tdata[i]['bott_den'][0::tran_space_interval,:].shape[0]):
    tmp[j,:] = mc.get_time_avg_data(da_mag[i][j,:]/sqrt(9.8*(bott_den[i][j,:]-surf_den[i][j,:])/1024*dpth[j,:]), dt, interval)
  tr_froude.append(tmp)

# output transect data:

for i in range(len(transects)):
  meccs_file_util.write_raw_data(out_base, "transect_"+tnames[i]+"_salt_intrusion_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 salt_intrusion[i].T, ["salinity intrusion"], "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "meters")
  meccs_file_util.write_raw_data(out_base, "transect_"+tnames[i]+"_output_xy_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 out_xy[i], ["x", "y"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "N/A", "meters")
  meccs_file_util.write_raw_data(out_base, "transect_"+tnames[i]+"_surf_vel_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 surf_vel[i].T, ["surface velocity"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "m/s")
  meccs_file_util.write_raw_data(out_base, "transect_"+tnames[i]+"_ac_davg_vel_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 da_vel[i].T, ["along channel depth averaged velocity"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "m/s")
  meccs_file_util.write_raw_data(out_base, "transect_"+tnames[i]+"_mag_depth_avg_vel_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 da_mag[i].T, ["magnitude of depth averaged velocity"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "m/s")
  meccs_file_util.write_raw_data(out_base, "transect_"+tnames[i]+"_surface_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 surf_salt[i].T, ["surface salinity"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "psu")
  meccs_file_util.write_raw_data(out_base, "transect_"+tnames[i]+"_bottom_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 bott_salt[i].T, ["bottom salinity"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "psu")
  meccs_file_util.write_raw_data(out_base, "transect_"+tnames[i]+"_davg_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 davg_salt[i].T, ["depth averaged salinity"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "psu")
  meccs_file_util.write_raw_data(out_base, "transect_"+tnames[i]+"_surface_rho_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 surf_den[i].T, ["surface density"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "kg/m^3")
  meccs_file_util.write_raw_data(out_base, "transect_"+tnames[i]+"_bottom_rho_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 bott_den[i].T, ["bottom density"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "kg/m^3")
  meccs_file_util.write_raw_data(out_base, "transect_"+tnames[i]+"_bottom_rho_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 bott_den[i].T, ["bottom density"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "kg/m^3")
  meccs_file_util.write_raw_data(out_base, "transect_"+tnames[i]+"_davg_rho_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 davg_den[i].T, ["depth averaged density"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "kg/m^3")
  meccs_file_util.write_raw_data(out_base, "transect_"+tnames[i]+"_froude_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tr_froude[i].T, ["interfacial Froude #"], "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "kg/m^3")


#tidal_amplitude
p_amps = []
for i in range(len(points)):
  p_amps.append(mc.get_amplitude(ptdata['elev'][i,:], dt, interval))

tmp = zeros((p_amps[0].shape[0], len(p_amps)), float)
for i in range(len(points)):
  tmp[:,i] = p_amps[i]

meccs_file_util.write_raw_data(out_base, "points_tidal_amplitude_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "meters")                  


#mean_river_flux
mean_flux = mc.get_time_avg_data(rvr['river_flow'][0,:], dt, interval)
meccs_file_util.write_raw_data(out_base, "mean_river_flux_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 mean_flux, ["mean river flux"], "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "m^3/s")                  

#mean volume calculations
mean_vol = []
prism = []
v_int_salt = []
v_avg_salt = []

for i in range(vdata['volume'].shape[0]):
  mean_vol.append(mc.get_time_avg_data(vdata['volume'][i,:], dt, interval))
  prism.append(mc.get_amplitude(vdata['volume'][i,:], dt, interval))
  
  v_int_salt.append(mc.get_time_avg_data(vdata['vol_int_salt'][i,:], dt, interval))
  v_avg_salt.append(mc.get_time_avg_data(vdata['vol_avg_salt'][i,:], dt, interval))

tmp = zeros((mean_vol[0].shape[0], len(regions)), float)
for i in range(len(regions)):
  tmp[:,i] = mean_vol[i]
meccs_file_util.write_raw_data(out_base, "volume_mean_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, rnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "m^3")

tmp = zeros((prism[0].shape[0], len(regions)), float)
for i in range(len(regions)):
  tmp[:,i] = prism[i]
meccs_file_util.write_raw_data(out_base, "volume_tidal_prism_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, rnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "m^3")

tmp = zeros((v_int_salt[0].shape[0], len(regions)), float)
for i in range(len(regions)):
  tmp[:,i] = v_int_salt[i]
meccs_file_util.write_raw_data(out_base, "volume_integrated_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, rnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu * m^3")

#volume calculations

flow_ratio = []
tidal_ex = []
fill_time = []
fw_flush = []
fsw = []

for i in range(vdata['volume'].shape[0]):
  #flow ratio: mean_river_flux * tidal_period / tidal_prism
  flow_ratio.append((mean_flux*tidal_period) / prism[i])
  #tidal exchange: msl_volume / tidal_prism
  tidal_ex.append(mean(mean_vol[i]) / prism[i])  
  #filling time: msl_volume / (mean_river_flux * tidal_period)
  fill_time.append(mean(mean_vol[i]) / (mean_flux*interval))
  #freshwater flushing time: volume_freshwater / mean_river_flux
  fw_flush.append((v_avg_salt[i]/ocean_salt)*mean_vol[i] / (mean_flux*interval))
  #fraction of saltwater
  fsw.append((v_avg_salt[i]/ocean_salt))

tmp = zeros((flow_ratio[0].shape[0], len(regions)), float)
for i in range(len(regions)):
  tmp[:,i] = flow_ratio[i]
meccs_file_util.write_raw_data(out_base, "volume_flow_ratio_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, rnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "dimensionless")

tmp = zeros((tidal_ex[0].shape[0], len(regions)), float)
for i in range(len(regions)):
  tmp[:,i] = tidal_ex[i]
meccs_file_util.write_raw_data(out_base, "volume_tidal_exchange_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, rnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "tidal periods")

tmp = zeros((fill_time[0].shape[0], len(regions)), float)
for i in range(len(regions)):
  tmp[:,i] = fill_time[i]
meccs_file_util.write_raw_data(out_base, "volume_filling_time_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, rnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "tidal periods")

tmp = zeros((fw_flush[0].shape[0], len(regions)), float)
for i in range(len(regions)):
  tmp[:,i] = fw_flush[i]
meccs_file_util.write_raw_data(out_base, "volume_fresh_water_flushing_time_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, rnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "seconds")

tmp = zeros((fsw[0].shape[0], len(regions)), float)
for i in range(len(regions)):
  tmp[:,i] = fsw[i]
meccs_file_util.write_raw_data(out_base, "volume_fraction_salt_water_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, rnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "% salt water")

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



#point calculations

HnR = []
froude = []
ac_hvel = []
elev = []
salt = []
s_salt = []
b_salt = []
ac_surf_vel = []
temp = []
b_temp = []
s_temp = []
rho = []
s_rho = []
b_rho = []


for i in range(ptdata['davg_salt'].shape[0]):
  HnR.append(mc.get_Hansen_Rattray(ptdata['surf_salt'][i,:], ptdata['bott_salt'][i,:], ptdata['davg_salt'][i,:], \
                                   ptdata['surf_hvel_dir'][i,:], rvr['river_flow'][0,:]/ptdata['depths'][0,i], dt, interval))

  froude.append(mc.get_time_avg_data(ptdata['davg_hvel_mag'][i,:]/sqrt(9.8*(ptdata['bott_den'][i,:]-ptdata['surf_den'][i,:])/1024*ptdata['depths'][0,i]), \
                                     dt, interval))

  ac_hvel.append(ptdata['davg_hvel_dir'][i,:])
  elev.append(ptdata['elev'][i,:])
  salt.append(ptdata['davg_salt'][i,:])
  s_salt.append(ptdata['surf_salt'][i,:])
  b_salt.append(ptdata['bott_salt'][i,:])
  ac_surf_vel.append(ptdata['surf_hvel_dir'][i,:])
  temp.append(ptdata['davg_temp'][i,:])
  s_temp.append(ptdata['surf_temp'][i,:])
  b_temp.append(ptdata['bott_temp'][i,:])
  rho.append(ptdata['davg_den'][i,:])
  s_rho.append(ptdata['surf_den'][i,:])
  b_rho.append(ptdata['bott_den'][i,:])  

tmp = zeros((HnR[0].shape[1], len(pnames)*2), float)
tmpnames = []
for i in range(len(pnames)):
  tmp[:,i*2:i*2+2] = HnR[i].T
  tmpnames.append("["+pnames[i]+": x\t y]")
meccs_file_util.write_raw_data(out_base, "point_Hansen_n_Rattray_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, tmpnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "dimensionless")

tmp = zeros((froude[0].shape[0], len(pnames)), float)
for i in range(len(pnames)):
  tmp[:,i] = froude[i]
meccs_file_util.write_raw_data(out_base, "point_froude_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "dimensionless")

tmp = zeros((ac_hvel[0].shape[0], len(pnames)), float)
for i in range(len(pnames)):
  tmp[:,i] = ac_hvel[i]
meccs_file_util.write_raw_data(out_base, "point_depth_avg_channel_velocity_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "m/s")

tmp = zeros((ac_hvel[0].shape[0], len(pnames)), float)
for i in range(len(pnames)):
  tmp[:,i] = ac_surf_vel[i]
meccs_file_util.write_raw_data(out_base, "point_surf_channel_velocity_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "m/s")

tmp = zeros((elev[0].shape[0], len(pnames)), float)
for i in range(len(pnames)):
  tmp[:,i] = elev[i]
meccs_file_util.write_raw_data(out_base, "point_elevation_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "m from MSL")

tmp = zeros((salt[0].shape[0], len(pnames)), float)
for i in range(len(pnames)):
  tmp[:,i] = salt[i]
meccs_file_util.write_raw_data(out_base, "point_depth_avg_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((salt[0].shape[0], len(pnames)), float)
for i in range(len(pnames)):
  tmp[:,i] = s_salt[i]
meccs_file_util.write_raw_data(out_base, "point_surf_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((salt[0].shape[0], len(pnames)), float)
for i in range(len(pnames)):
  tmp[:,i] = b_salt[i]
meccs_file_util.write_raw_data(out_base, "point_bott_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((temp[0].shape[0], len(pnames)), float)
for i in range(len(pnames)):
  tmp[:,i] = temp[i]
meccs_file_util.write_raw_data(out_base, "point_depth_avg_temp_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((temp[0].shape[0], len(pnames)), float)
for i in range(len(pnames)):
  tmp[:,i] = s_temp[i]
meccs_file_util.write_raw_data(out_base, "point_surf_temp_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((temp[0].shape[0], len(pnames)), float)
for i in range(len(pnames)):
  tmp[:,i] = b_temp[i]
meccs_file_util.write_raw_data(out_base, "point_bott_temp_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((salt[0].shape[0], len(pnames)), float)
for i in range(len(pnames)):
  tmp[:,i] = rho[i]
meccs_file_util.write_raw_data(out_base, "point_depth_avg_den_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((salt[0].shape[0], len(pnames)), float)
for i in range(len(pnames)):
  tmp[:,i] = s_rho[i]
meccs_file_util.write_raw_data(out_base, "point_surf_den_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((salt[0].shape[0], len(pnames)), float)
for i in range(len(pnames)):
  tmp[:,i] = b_rho[i]
meccs_file_util.write_raw_data(out_base, "point_bott_den_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

