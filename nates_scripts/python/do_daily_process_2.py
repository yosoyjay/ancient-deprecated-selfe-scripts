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

#data types for individual point values
#data_types = ["salt", "temp", "hvel_u", "hvel_v", "elev"]
data_types = ["salt", "temp", "hvel_u", "hvel_v", "elev", "zcor"]

#for computing significant equality
small = 0.0001

# set the time (which also dictates output file name)
start_time = ydhms2corie(yr, day, 0, 15, 0)
end_time = start_time + ndays-0.25/24  #puts it right at the end of the day

sp = shape_procs.shape_procs()

# read all the shapes in: assume they are in the shape dir and named "points.txt", "cross_sections.txt", "transects.txt" and "regions.txt"
if os.access(shape_dir+"/points.txt", os.F_OK):
  (points, pnames, privers) = meccs_file_util.read_cs_file(shape_dir+"/points.txt", 1)   #uses cs file format
  do_points=1
else:
  do_points=0
if os.access(shape_dir+"/cross_sections.txt", os.F_OK):  
  (cross_sections, cnames, crivers) = meccs_file_util.read_cs_file(shape_dir+"/cross_sections.txt", 1)
  do_cs=1
else:
  do_cs=0  
if os.access(shape_dir+"/transects.txt", os.F_OK):
#if 0:
  (transects, tnames, trivers) = meccs_file_util.read_cs_file(shape_dir+"/transects.txt", 1)
  do_transects=1
else:
  do_transects=0  
if os.access(shape_dir+"/regions.txt", os.F_OK):  
  (regions, rnames, rrivers) = meccs_file_util.read_cs_file(shape_dir+"/regions.txt", 0)
  do_regions=1
else:
  do_regions=0

#pdb.set_trace()

for i in range(ndays):
  curr_time = start_time+i
  in_dir = "%s/%s/%d-%03d/process/" % (base_dir, estuaryid, yr, day+i)
  if (curr_time==start_time):
    rvr = sp.read_cs_data(in_dir+"/river_flow.dat")
  else:
    tmp = sp.read_cs_data(in_dir+"/river_flow.dat")
    rvr['river_flow'] = concatenate((rvr['river_flow'], tmp['river_flow']), 0)
    rvr['end_time'] = tmp['end_time']

if (do_transects):
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

if (do_cs):
  for i in range(ndays):  
    curr_time = start_time+i
    in_dir = "%s/%s/%d-%03d/process/" % (base_dir, estuaryid, yr, day+i)        
    if (i==0):
      csdata = sp.read_cs_data(in_dir+"/cross_section_data.dat")
    else:
      tmp=sp.read_cs_data(in_dir+"/cross_section_data.dat")
      for akey in tmp.keys():
        if (akey!='dt' and akey!='name' and akey!='end_time' and akey!='start_time'):        
          csdata[akey] = concatenate((csdata[akey], tmp[akey]), 1)
      csdata['end_time'] = tmp['end_time']

if (do_points):
  ipts = []
  for i in range(ndays):  
    curr_time = start_time+i
    in_dir = "%s/%s/%d-%03d/process/" % (base_dir, estuaryid, yr, day+i)        
    if (i==0):
      ptdata = sp.read_cs_data(in_dir+"/point_data.dat")
    else:
      tmp = sp.read_cs_data(in_dir+"/point_data.dat")
      for akey in tmp.keys():
        if (akey!='dt' and akey!='names' and akey!='end_time' and akey!='start_time' and akey!='depths'):        
          ptdata[akey] = concatenate((ptdata[akey], tmp[akey]), 1)
      ptdata['end_time'] = tmp['end_time']

    for j in range(len(pnames)):
      if (i==0):
        ipts.append({})
#        ipts[j] = {}
        ipts[j]["zmsl"] = meccs_file_util.read_raw_data(in_dir+pnames[j]+"_zmsl.dat")
      for dtype in data_types:
        fname = in_dir+pnames[j]+"_"+dtype+".dat"
        if (i==0):
          ipts[j][dtype] = meccs_file_util.read_raw_data(fname)
        else:
          ipts[j][dtype] = meccs_file_util.col_stack_data([ipts[j][dtype], meccs_file_util.read_raw_data(fname)])
  try:
    for j in range(len(pnames)):
      if (ipts[j]['elev']['data'].ndim > 1):
        for j in range(len(pnames)):
          ipts[j]['elev']['data'] = ipts[j]['elev']['data'].T.reshape(ipts[j]['elev']['data'].shape[0]*ipts[j]['elev']['data'].shape[1])
  except AttributeError:
    pass
  
if (do_regions):
  for i in range(ndays):  
    curr_time = start_time+i
    in_dir = "%s/%s/%d-%03d/process/" % (base_dir, estuaryid, yr, day+i)        
    if (i==0):
      vdata = sp.read_cs_data(in_dir+"/volume_data.dat")
    else:
      tmp = sp.read_cs_data(in_dir+"/volume_data.dat")
      for akey in tmp.keys():
        if (akey!='dt' and akey!='name' and akey!='end_time' and akey!='start_time'):        
          vdata[akey] = concatenate((vdata[akey], tmp[akey]), 1)
      vdata['end_time'] = tmp['end_time']

interval = 24.8 * 3600 *3
tidal_period = 12.4 * 3600
dt = 900
mc = meccs.meccs()
ocean_salt = 33.5
ocean_den = 1024.0
salt_intr = 2.0 


##transect data - do every 10 pts
#tran_space_interval = 10
#salt_intrusion = []
#tr_froude = []
#surf_vel = []
#da_vel = []
#da_mag = []
#surf_salt = []
#bott_salt = []
#surf_den = []
#bott_den = []
#davg_den = []
#davg_salt = []
#out_xy = []
#
#if (do_transects):
#  for i in range(len(transects)):
#    sp.add_shape(transects[i], "transect", 50, tnames[i])
#    temp_shape = sp.get_shape(tnames[i])
#  
#    salt_intrusion.append(mc.get_salinity_intrusion(temp_shape, tdata[i]['dmax_salt'], dt, interval, salt_intr))
#    out_xy.append(temp_shape[0::tran_space_interval, 0:2])
#    surf_vel.append(tdata[i]['surf_hvel_dir'][0::tran_space_interval,:])
#    da_vel.append(tdata[i]['davg_hvel_dir'][0::tran_space_interval,:])
#    da_mag.append(tdata[i]['davg_hvel_mag'][0::tran_space_interval,:])
#    surf_salt.append(tdata[i]['surf_salt'][0::tran_space_interval,:])
#    bott_salt.append(tdata[i]['bott_salt'][0::tran_space_interval,:])
#    davg_salt.append(tdata[i]['davg_salt'][0::tran_space_interval,:])
#    surf_den.append(tdata[i]['surf_den'][0::tran_space_interval,:])
#    bott_den.append(tdata[i]['bott_den'][0::tran_space_interval,:])  
#    davg_den.append(tdata[i]['davg_den'][0::tran_space_interval,:])
#  
#    tstep_int = ceil(interval/dt)  #keeping it simple and slightly innacurate
#    tmp = zeros((tdata[i]['bott_den'][0::tran_space_interval,:].shape[0], int(tdata[i]['bott_den'].shape[1]/tstep_int)), float)
#    dpth = tdata[i]['depth'][0::tran_space_interval,:]
#    for j in range(tdata[i]['bott_den'][0::tran_space_interval,:].shape[0]):
#      tmp[j,:] = mc.get_time_avg_data(da_mag[i][j,:]/sqrt(9.8*(bott_den[i][j,:]-surf_den[i][j,:])/1024*dpth[j,:]), dt, interval)
#    tr_froude.append(tmp)
#  
#  # output transect data:
#  
#  for i in range(len(transects)):
#    meccs_file_util.write_raw_data(out_dir, "transect_"+tnames[i]+"_salt_intrusion_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
#                                   salt_intrusion[i].T, ["salinity intrusion"], "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "meters")
#    meccs_file_util.write_raw_data(out_dir, "transect_"+tnames[i]+"_output_xy_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
#                                   out_xy[i], ["x", "y"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "N/A", "meters")
#    meccs_file_util.write_raw_data(out_dir, "transect_"+tnames[i]+"_surf_vel_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
#                                   surf_vel[i].T, ["surface velocity"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "m/s")
#    meccs_file_util.write_raw_data(out_dir, "transect_"+tnames[i]+"_ac_davg_vel_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
#                                   da_vel[i].T, ["along channel depth averaged velocity"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "m/s")
#    meccs_file_util.write_raw_data(out_dir, "transect_"+tnames[i]+"_mag_depth_avg_vel_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
#                                   da_mag[i].T, ["magnitude of depth averaged velocity"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "m/s")
#    meccs_file_util.write_raw_data(out_dir, "transect_"+tnames[i]+"_surface_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
#                                   surf_salt[i].T, ["surface salinity"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "psu")
#    meccs_file_util.write_raw_data(out_dir, "transect_"+tnames[i]+"_bottom_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
#                                   bott_salt[i].T, ["bottom salinity"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "psu")
#    meccs_file_util.write_raw_data(out_dir, "transect_"+tnames[i]+"_davg_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
#                                   davg_salt[i].T, ["depth averaged salinity"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "psu")
#    meccs_file_util.write_raw_data(out_dir, "transect_"+tnames[i]+"_surface_rho_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
#                                   surf_den[i].T, ["surface density"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "kg/m^3")
#    meccs_file_util.write_raw_data(out_dir, "transect_"+tnames[i]+"_bottom_rho_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
#                                   bott_den[i].T, ["bottom density"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "kg/m^3")
#    meccs_file_util.write_raw_data(out_dir, "transect_"+tnames[i]+"_bottom_rho_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
#                                   bott_den[i].T, ["bottom density"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "kg/m^3")
#    meccs_file_util.write_raw_data(out_dir, "transect_"+tnames[i]+"_davg_rho_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
#                                   davg_den[i].T, ["depth averaged density"], "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "kg/m^3")
#    meccs_file_util.write_raw_data(out_dir, "transect_"+tnames[i]+"_froude_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
#                                   tr_froude[i].T, ["interfacial Froude #"], "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "kg/m^3")


##tidal_amplitude
#p_amps = []
#for i in range(len(points)):
#  p_amps.append(mc.get_amplitude(ptdata['elev'][i,:], dt, interval))
#
#tmp = zeros((p_amps[0].shape[0], len(p_amps)), float)
#for i in range(len(points)):
#  tmp[:,i] = p_amps[i]
#
#meccs_file_util.write_raw_data(out_dir, "points_tidal_amplitude_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
#                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "meters")                  
#
#
##mean_river_flux
##pdb.set_trace()
#
#if (rvr['river_flow'].ndim > 1):
#  mean_flux = zeros((floor(rvr['river_flow'].shape[0]*dt/interval), rvr['river_flow'].shape[1]), float)
#  for i in range(rvr['river_flow'].shape[1]):
#    mean_flux[:,i] = mc.get_time_avg_data(rvr['river_flow'][:,i], dt, interval)
#else:
#  mean_flux = mc.get_time_avg_data(rvr['river_flow'], dt, interval)
#
##mean_flux = mc.get_time_avg_data(rvr['river_flow'][0,:], dt, interval)
#meccs_file_util.write_raw_data(out_dir, "mean_river_flux_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
#                               mean_flux, ["mean river flux"], "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "m^3/s")                  
#
##mean volume calculations
#mean_vol = []
#prism = []
#v_int_salt = []
#v_avg_salt = []
#
#for i in range(vdata['volume'].shape[0]):
#  mean_vol.append(mc.get_time_avg_data(vdata['volume'][i,:], dt, interval))
#  prism.append(mc.get_amplitude(vdata['volume'][i,:], dt, interval))  
#  v_int_salt.append(mc.get_time_avg_data(vdata['vol_int_salt'][i,:], dt, interval))
#  v_avg_salt.append(mc.get_time_avg_data(vdata['vol_avg_salt'][i,:], dt, interval))
#
#tmp = zeros((mean_vol[0].shape[0], len(regions)), float)
#for i in range(len(regions)):
#  tmp[:,i] = mean_vol[i]
#meccs_file_util.write_raw_data(out_dir, "volume_mean_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
#                                 tmp, rnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "m^3")
#
#tmp = zeros((prism[0].shape[0], len(regions)), float)
#for i in range(len(regions)):
#  tmp[:,i] = prism[i]
#meccs_file_util.write_raw_data(out_dir, "volume_tidal_prism_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
#                                 tmp, rnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "m^3")
#
#tmp = zeros((v_int_salt[0].shape[0], len(regions)), float)
#for i in range(len(regions)):
#  tmp[:,i] = v_int_salt[i]
#meccs_file_util.write_raw_data(out_dir, "volume_integrated_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
#                                 tmp, rnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu * m^3")
#
##volume calculations
#
#flow_ratio = []
#tidal_ex = []
#fill_time = []
#fw_flush = []
#fsw = []
#
#pdb.set_trace()
#
#for i in range(vdata['volume'].shape[0]):
#  mflux = sum(mean_flux[:,rrivers[i]],1)
#  #flow ratio: mean_river_flux * tidal_period / tidal_prism
#  flow_ratio.append((mflux*tidal_period) / prism[i])
#  #tidal exchange: msl_volume / tidal_prism
#  tidal_ex.append(mean(mean_vol[i]) / prism[i])  
#  #filling time: msl_volume / (mean_river_flux * tidal_period)
#  fill_time.append(mean(mean_vol[i]) / (mflux*interval))
#  #freshwater flushing time: volume_freshwater / mean_river_flux
#  fw_flush.append((v_avg_salt[i]/ocean_salt)*mean_vol[i] / (mflux*interval))
#  #fraction of saltwater
#  fsw.append((v_avg_salt[i]/ocean_salt))
#
#tmp = zeros((flow_ratio[0].shape[0], len(regions)), float)
#for i in range(len(regions)):
#  tmp[:,i] = flow_ratio[i]
#meccs_file_util.write_raw_data(out_dir, "volume_flow_ratio_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
#                                 tmp, rnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "dimensionless")
#
#tmp = zeros((tidal_ex[0].shape[0], len(regions)), float)
#for i in range(len(regions)):
#  tmp[:,i] = tidal_ex[i]
#meccs_file_util.write_raw_data(out_dir, "volume_tidal_exchange_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
#                                 tmp, rnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "tidal periods")
#
#tmp = zeros((fill_time[0].shape[0], len(regions)), float)
#for i in range(len(regions)):
#  tmp[:,i] = fill_time[i]
#meccs_file_util.write_raw_data(out_dir, "volume_filling_time_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
#                                 tmp, rnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "tidal periods")
#
#tmp = zeros((fw_flush[0].shape[0], len(regions)), float)
#for i in range(len(regions)):
#  tmp[:,i] = fw_flush[i]
#meccs_file_util.write_raw_data(out_dir, "volume_fresh_water_flushing_time_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
#                                 tmp, rnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "seconds")
#
#tmp = zeros((fsw[0].shape[0], len(regions)), float)
#for i in range(len(regions)):
#  tmp[:,i] = fsw[i]
#meccs_file_util.write_raw_data(out_dir, "volume_fraction_salt_water_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
#                                 tmp, rnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "% salt water")

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

#pdb.set_trace()

for i in range(csdata['davg_salt'].shape[0]):
  rflow = sum(rvr['river_flow'][:,crivers[i]], 1)
  
#  cs_HnR.append(mc.get_Hansen_Rattray(csdata['cs_surf_salt'][i,:], csdata['cs_bott_salt'][i,:], csdata['davg_salt'][i,:], \
#                                      csdata['cs_surf_hvel_dir'][i,:], rvr['river_flow'][0,:]/csdata['cs_areas'][0,:], dt, interval))

#  temp = mc.get_Hansen_Rattray(csdata['cs_surf_salt'][i,:], csdata['cs_bott_salt'][i,:], csdata['davg_salt'][i,:], csdata['cs_surf_hvel_dir'][i,:], rflow/csdata['cs_areas'][0,:], dt, interval)

#  cs_HnR.append(temp);

  cs_HnR.append(mc.get_Hansen_Rattray(csdata['cs_surf_salt'][i,:], csdata['cs_bott_salt'][i,:], csdata['davg_salt'][i,:], csdata['cs_surf_hvel_dir'][i,:], rflow/csdata['cs_areas'][0,:], dt, interval))

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
meccs_file_util.write_raw_data(out_dir, "csection_Hansen_n_Rattray_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, tmpnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "dimensionless")



tmp = zeros((cs_froude[0].shape[0], len(cnames)), float)
for i in range(len(cnames)):
  tmp[:,i] = cs_froude[i]
meccs_file_util.write_raw_data(out_dir, "csection_froude_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, cnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "dimensionless")

tmp = zeros((cs_flux[0].shape[0], len(cnames)), float)
for i in range(len(cnames)):
  tmp[:,i] = cs_flux[i]
meccs_file_util.write_raw_data(out_dir, "csection_flux_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, cnames, "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "m^3/s")

tmp = zeros((cs_area[0].shape[0], len(cnames)), float)
for i in range(len(cnames)):
  tmp[:,i] = cs_area[i]
meccs_file_util.write_raw_data(out_dir, "csection_area_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, cnames, "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "m^2")

tmp = zeros((cs_salt[0].shape[0], len(cnames)), float)
for i in range(len(cnames)):
  tmp[:,i] = cs_salt[i]
meccs_file_util.write_raw_data(out_dir, "csection_avg_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, cnames, "day:\t"+str(day)+"\tyear\t"+str(yr), dt, "seconds", "psu")

tmp = zeros((cs_salt[0].shape[0], len(cnames)), float)
for i in range(len(cnames)):
  tmp[:,i] = cs_s_salt[i]
meccs_file_util.write_raw_data(out_dir, "csection_surf_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, cnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((cs_salt[0].shape[0], len(cnames)), float)
for i in range(len(cnames)):
  tmp[:,i] = cs_b_salt[i]
meccs_file_util.write_raw_data(out_dir, "csection_bott_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, cnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((cs_temp[0].shape[0], len(cnames)), float)
for i in range(len(cnames)):
  tmp[:,i] = cs_temp[i]
meccs_file_util.write_raw_data(out_dir, "csection_depth_avg_temp_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, cnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((cs_temp[0].shape[0], len(cnames)), float)
for i in range(len(cnames)):
  tmp[:,i] = cs_s_temp[i]
meccs_file_util.write_raw_data(out_dir, "csection_surf_temp_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, cnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((cs_temp[0].shape[0], len(cnames)), float)
for i in range(len(cnames)):
  tmp[:,i] = cs_b_temp[i]
meccs_file_util.write_raw_data(out_dir, "csection_bott_temp_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, cnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((cs_salt[0].shape[0], len(cnames)), float)
for i in range(len(cnames)):
  tmp[:,i] = cs_rho[i]
meccs_file_util.write_raw_data(out_dir, "csection_depth_avg_den_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, cnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((cs_salt[0].shape[0], len(cnames)), float)
for i in range(len(cnames)):
  tmp[:,i] = cs_s_rho[i]
meccs_file_util.write_raw_data(out_dir, "csection_surf_den_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, cnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((cs_salt[0].shape[0], len(cnames)), float)
for i in range(len(cnames)):
  tmp[:,i] = cs_b_rho[i]
meccs_file_util.write_raw_data(out_dir, "csection_bott_den_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
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
#  HnR.append(mc.get_Hansen_Rattray(ptdata['surf_salt'][i,:], ptdata['bott_salt'][i,:], ptdata['davg_salt'][i,:], \
#                                   ptdata['surf_hvel_dir'][i,:], rvr['river_flow'][0,:]/ptdata['depths'][0,i], dt, interval))

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

#tmp = zeros((HnR[0].shape[1], len(pnames)*2), float)
#tmpnames = []
#for i in range(len(pnames)):
#  tmp[:,i*2:i*2+2] = HnR[i].T
#  tmpnames.append("["+pnames[i]+": x\t y]")
#meccs_file_util.write_raw_data(out_dir, "point_Hansen_n_Rattray_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
#                                 tmp, tmpnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "dimensionless")

tmp = zeros((froude[0].shape[0], len(pnames)), float)
for i in range(len(pnames)):
  tmp[:,i] = froude[i]
meccs_file_util.write_raw_data(out_dir, "point_froude_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "dimensionless")

tmp = zeros((ac_hvel[0].shape[0], len(pnames)), float)
for i in range(len(pnames)):
  tmp[:,i] = ac_hvel[i]
meccs_file_util.write_raw_data(out_dir, "point_depth_avg_channel_velocity_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "m/s")

tmp = zeros((ac_hvel[0].shape[0], len(pnames)), float)
for i in range(len(pnames)):
  tmp[:,i] = ac_surf_vel[i]
meccs_file_util.write_raw_data(out_dir, "point_surf_channel_velocity_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "m/s")

tmp = zeros((elev[0].shape[0], len(pnames)), float)
for i in range(len(pnames)):
  tmp[:,i] = elev[i]
meccs_file_util.write_raw_data(out_dir, "point_elevation_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "m from MSL")

tmp = zeros((salt[0].shape[0], len(pnames)), float)
for i in range(len(pnames)):
  tmp[:,i] = salt[i]
meccs_file_util.write_raw_data(out_dir, "point_depth_avg_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((salt[0].shape[0], len(pnames)), float)
for i in range(len(pnames)):
  tmp[:,i] = s_salt[i]
meccs_file_util.write_raw_data(out_dir, "point_surf_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((salt[0].shape[0], len(pnames)), float)
for i in range(len(pnames)):
  tmp[:,i] = b_salt[i]
meccs_file_util.write_raw_data(out_dir, "point_bott_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((temp[0].shape[0], len(pnames)), float)
for i in range(len(pnames)):
  tmp[:,i] = temp[i]
meccs_file_util.write_raw_data(out_dir, "point_depth_avg_temp_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((temp[0].shape[0], len(pnames)), float)
for i in range(len(pnames)):
  tmp[:,i] = s_temp[i]
meccs_file_util.write_raw_data(out_dir, "point_surf_temp_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((temp[0].shape[0], len(pnames)), float)
for i in range(len(pnames)):
  tmp[:,i] = b_temp[i]
meccs_file_util.write_raw_data(out_dir, "point_bott_temp_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((salt[0].shape[0], len(pnames)), float)
for i in range(len(pnames)):
  tmp[:,i] = rho[i]
meccs_file_util.write_raw_data(out_dir, "point_depth_avg_den_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((salt[0].shape[0], len(pnames)), float)
for i in range(len(pnames)):
  tmp[:,i] = s_rho[i]
meccs_file_util.write_raw_data(out_dir, "point_surf_den_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")

tmp = zeros((salt[0].shape[0], len(pnames)), float)
for i in range(len(pnames)):
  tmp[:,i] = b_rho[i]
meccs_file_util.write_raw_data(out_dir, "point_bott_den_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 tmp, pnames, "day:\t"+str(day)+"\tyear\t"+str(yr), interval, "seconds", "psu")
  


############## do data for individual points ####################

#for absolute depth interpolation
out_depth_abs = array([1.0])
abs_od_str = ""
for i in range(out_depth_abs.shape[0]):
  abs_od_str = abs_od_str+" "+str(out_depth_abs[i])

#for ratio depth interpolation
out_depth_ratios = array([1, 0.75, 0.5, 0.25, 0])
od_str = ""
for i in range(out_depth_ratios.shape[0]):
  od_str = od_str+" "+str(out_depth_ratios[i])

for i in range(len(ipts)):
  ret_salt_a = interpolate.interp_z_abs(ipts[i]['salt']['data'], out_depth_abs, ipts[i]['zcor']['data'])
  ret_temp_a = interpolate.interp_z_abs(ipts[i]['temp']['data'], out_depth_abs, ipts[i]['zcor']['data'])
  ret_hvel_a = interpolate.interp_z_abs(ipts[i]['hvel_u']['data'], out_depth_abs, ipts[i]['zcor']['data'])\
                     * points[i][0,2] + \
                     interpolate.interp_z_abs(ipts[i]['hvel_v']['data'], out_depth_abs, ipts[i]['zcor']['data'])\
                     * points[i][0,3]

  zmsl = ipts[i]['zmsl']['data']
  ret_salt = interpolate.interp_z_ratio(ipts[i]['salt']['data'], out_depth_ratios, -abs(ptdata['depths'][0,i]), zmsl[:,0])
  ret_temp = interpolate.interp_z_ratio(ipts[i]['temp']['data'], out_depth_ratios, -abs(ptdata['depths'][0,i]), zmsl[:,0])
  ret_hvel = interpolate.interp_z_ratio(ipts[i]['hvel_u']['data'], out_depth_ratios, -abs(ptdata['depths'][0,i]), zmsl[:,0])\
                     * points[i][0,2] + \
                     interpolate.interp_z_ratio(ipts[i]['hvel_v']['data'], out_depth_ratios, -abs(ptdata['depths'][0,i]), zmsl[:,0])\
                     * points[i][0,3]
  ret_elev = ipts[i]['elev']['data']
  ret_elev = ret_elev.reshape(1, ret_elev.shape[0])

  description = od_str + " of depth " + str(ptdata['depths'][0,i]) + " from MSL"
  meccs_file_util.write_raw_data(out_dir, pnames[i]+"_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 ret_salt, description, "day:\t"+str(day)+"\tyear\t"+str(yr), "900", "seconds", "psu")
  meccs_file_util.write_csv(out_dir+pnames[i]+"_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".csv", ret_salt)

  meccs_file_util.write_raw_data(out_dir, pnames[i]+"_temp_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 ret_temp, description, "day:\t"+str(day)+"\tyear\t"+str(yr), "900", "seconds", "psu")
  meccs_file_util.write_csv(out_dir+pnames[i]+"_temp_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".csv", ret_temp)

  meccs_file_util.write_raw_data(out_dir, pnames[i]+"_hvel_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 ret_hvel, description, "day:\t"+str(day)+"\tyear\t"+str(yr), "900", "seconds", "psu")
  meccs_file_util.write_csv(out_dir+pnames[i]+"_hvel_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".csv", ret_hvel)

  meccs_file_util.write_raw_data(out_dir, pnames[i]+"_elev_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 ret_elev, description, "day:\t"+str(day)+"\tyear\t"+str(yr), "900", "seconds", "psu")
  meccs_file_util.write_csv(out_dir+pnames[i]+"_elev_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".csv", ret_elev)

  meccs_file_util.write_raw_data(out_dir, pnames[i]+"_abs_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 ret_salt_a, description, "day:\t"+str(day)+"\tyear\t"+str(yr), "900", "seconds", "psu")
  meccs_file_util.write_csv(out_dir+pnames[i]+"_abs_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".csv", ret_salt_a)

  meccs_file_util.write_raw_data(out_dir, pnames[i]+"_abs_temp_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 ret_temp_a, description, "day:\t"+str(day)+"\tyear\t"+str(yr), "900", "seconds", "psu")
  meccs_file_util.write_csv(out_dir+pnames[i]+"_abs_temp_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".csv", ret_temp_a)

  meccs_file_util.write_raw_data(out_dir, pnames[i]+"_abs_hvel_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".dat", \
                                 ret_hvel_a, description, "day:\t"+str(day)+"\tyear\t"+str(yr), "900", "seconds", "psu")
  meccs_file_util.write_csv(out_dir+pnames[i]+"_abs_hvel_"+str(yr)+"_"+str(day)+"_"+str(ndays)+".csv", ret_hvel_a)

# metadata file

fid = open(out_dir+"raw_pt_metadata.txt", "w")
fid.write("point names:\n")
for i in range(len(pnames)):
  fid.write(pnames[i])
  if (i==len(pnames)-1):
    fid.write("\n")
  else:
    fid.write(", ")
fid.write("_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+"\n")  
fid.write("data type: salinity, depth type: ratios, depths: "+od_str+", start_time: "+str(start_time)+", interval: "+str(900./86400.)+"\n")
fid.write("_temp_"+str(yr)+"_"+str(day)+"_"+str(ndays)+"\n")  
fid.write("data type: temperature, depth type: ratios, depths: "+od_str+", start_time: "+str(start_time)+", interval: "+str(900./86400.)+"\n")
fid.write("_hvel_"+str(yr)+"_"+str(day)+"_"+str(ndays)+"\n")  
fid.write("data type: horizontal velocity, depth type: ratios, depths: "+od_str+", start_time: "+str(start_time)+", interval: "+str(900./86400.)+"\n")
fid.write("_elev_"+str(yr)+"_"+str(day)+"_"+str(ndays)+"\n")  
fid.write("data type: surface elevation, depth type: 2D, depths: 0, start_time: "+str(start_time)+", interval: "+str(900./86400.)+"\n")

fid.write("_abs_salt_"+str(yr)+"_"+str(day)+"_"+str(ndays)+"\n")  
fid.write("data type: salinity, depth type: absolute, depths: "+abs_od_str+", start_time: "+str(start_time)+", interval: "+str(900./86400.)+"\n")
fid.write("_abs_temp_"+str(yr)+"_"+str(day)+"_"+str(ndays)+"\n")  
fid.write("data type: temperature, depth type: absolute, depths: "+abs_od_str+", start_time: "+str(start_time)+", interval: "+str(900./86400.)+"\n")
fid.write("_abs_hvel_"+str(yr)+"_"+str(day)+"_"+str(ndays)+"\n")  
fid.write("data type: horizontal velocity, depth type: absolute, depths: "+abs_od_str+", start_time: "+str(start_time)+", interval: "+str(900./86400.)+"\n")

fid.close()
