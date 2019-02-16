import sys
sys.path.append('/home/users/hyde/ourestuaries/pythonlib/cmop/')
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
import os
import pdb

if len(sys.argv) < 7:
  sys.exit("usage: python do_daily_extract.py [estuary id] [base_dir] [shape_dir] [year] [day] [ndays] [full refresh]")

estuaryid  = sys.argv[1]
base_dir   = sys.argv[2]
shape_dir  = sys.argv[3]
yr         = int(sys.argv[4])
day        = int(sys.argv[5])
ndays      = int(sys.argv[6])
do_refresh = sys.argv[7]

# set the time (which also dictates output file name)
start_time = ydhms2corie(yr, day, 0, 15, 0)
end_time = start_time + ndays-0.25/24  #puts it right at the end of the day

#print "start_time: "+str(start_time)

mr = mreader(base_dir+"/"+estuaryid+"/", array([]), start_time, end_time, 0, 16, 1, 0)
sp = shape_procs.shape_procs()

## get the river flow data
#for i in range(ndays):
#  curr_time = start_time+i
#  river_flow=mr.get_river_flux(base_dir+"/"+estuaryid+"/", curr_time, curr_time+(1-0.25/24), 16, 900, 0)
#  data = {}
#  data["river_flow"] = river_flow
#  data["start_time"] = start_time
#  data["end_time"] = end_time
#  data["dt"] = 900
#  dname = "%s/%s/%d-%03d/process/" % (base_dir, estuaryid, yr, day+i)
#
#  if (not(os.path.exists(dname))):
#    os.mkdir(dname)
#  
#  sp.save_cs_data(data, dname+"/river_flow.dat")
#
## Cross sections data: flux, area, cross-sectionally averaged and normalized surface and bottom for temp, salt density and hvel (surf onl
#if (os.path.exists(shape_dir+"/cross_sections.txt")):
#  sp = shape_procs.shape_procs()
#
#  (cross_sections, cnames, crivers) = meccs_file_util.read_cs_file(shape_dir+"/cross_sections.txt", 1)
#  if len(cross_sections)>0:
#    for i in range(len(cross_sections)):
#      sp.add_shape(cross_sections[i], "cross_sec", 50, cnames[i])
#
#    if (meccs_file_util.file_exists(shape_dir+"/cross_sections.ob") and int(do_refresh)==0):
##      print "using cross section ob file "+shape_dir+"cross_sections.ob\n"      
#      sp.ob_file=shape_dir+"/cross_sections.ob"
#    else:
##      print "calculating new cross section ob file\n"
#      agr = Gr()
#      [year, day, ts] = model_util.corie2ydn(start_time)
#      dname = "%s/%d-%03d/run/" % (base_dir+"/"+estuaryid+"/", yr, day)      
#      agr.readHGrid(dname+"hgrid.gr3")
#      aob = Ob()
#      sppts = sp.get_pts("cross_sec")
#      
#      aob.grid_init(agr, array([sppts[0], sppts[1]]).T)
#      aob.save_ob(shape_dir+"/cross_sections.ob")
#      sp.ob_file = shape_dir+"/cross_sections.ob"
#
#    if (mean(abs(cross_sections[0][:,2]))==0 and mean(abs(cross_sections[0][:,3]))==0) or int(do_refresh)==1:
#      sp.gen_channel_norms(base_dir+"/"+estuaryid+"/", start_time, start_time+2, "cross_sec", 16, 1, 0)
#      meccs_file_util.write_cs_file(shape_dir+"/cross_sections.txt", sp.shapes, cnames, crivers)
#
#
#    for i in range(ndays):
#      curr_time = start_time+i
#
#      data = sp.gen_cross_section_data(base_dir+"/"+estuaryid+"/", curr_time, curr_time+(1-0.25/24), 16, "base", 1, 1, 0)
#      data["start_time"] = curr_time
#      data["end_time"] = curr_time+(1-0.25/24)
#      dname = "%s/%s/%d-%03d/process/" % (base_dir, estuaryid, yr, day+i)
#      sp.save_cs_data(data, dname+"/cross_section_data.dat")


# transect data    
if (os.path.exists(shape_dir+"/transects.txt")):
  (transects, tnames, trivers) = meccs_file_util.read_cs_file(shape_dir+"/transects.txt", 1)
  if len(transects)>0:
    for i in range(len(transects)):
      sp = shape_procs.shape_procs()      
      sp.add_shape(transects[i], "transect", 50, tnames[i])

      if (meccs_file_util.file_exists(shape_dir+"/"+tnames[i]+"_transects.ob") and int(do_refresh)==0):
        sp.ob_file=shape_dir+"/"+tnames[i]+"_transects.ob"
      else:
        agr = Gr()
        [year, day, ts] = model_util.corie2ydn(start_time)
        dname = "%s/%d-%03d/run/" % (base_dir+"/"+estuaryid+"/", yr, day)      
        agr.readHGrid(dname+"hgrid.gr3")
        aob = Ob()
        sppts = sp.get_pts("transect")
        
        aob.grid_init(agr, array([sppts[0], sppts[1]]).T)
        aob.save_ob(shape_dir+"/"+tnames[i]+"_transects.ob")
        sp.ob_file = shape_dir+"/"+tnames[i]+"_transects.ob"
  
      if (mean(abs(transects[0][:,2]))==0 and mean(abs(transects[0][:,3]))==0) or int(do_refresh)==1:
        sp.gen_channel_norms(base_dir+"/"+estuaryid+"/", start_time, start_time+2, "transect", 16, 1, 0)
#        meccs_file_util.write_cs_file(shape_dir+"/"+tnames[i]+"_transect.txt", sp.shapes, tnames, crivers)
#          have to work out how to store ac norms later: probably just move this chunk of code + appropriate constructors up
  
      for j in range(ndays):
        curr_time = start_time+j
  
        data = sp.gen_transect_data(base_dir+"/"+estuaryid+"/", curr_time, curr_time+(1-0.25/24), tnames[0], 16, 1, 1, 0)
        data["start_time"] = curr_time
        data["end_time"] = curr_time+(1-0.25/24)
        data["dt"] = 900
        dname = "%s/%s/%d-%03d/process/" % (base_dir, estuaryid, yr, day+j)
        sp.save_cs_data(data, dname+"/"+tnames[i]+"_transect_data.dat")

 #points
if (os.path.exists(shape_dir+"/points.txt")):
  sp = shape_procs.shape_procs()
  (points, pnames, privers) = meccs_file_util.read_cs_file(shape_dir+"/points.txt", 1)
  apoints = zeros((len(points), 4), float)
  if len(points)>0:
    for i in range(len(points)):
      sp.add_shape(points[i], "point", 0, pnames[i])
      apoints[i,:] = points[i][0]

    if (meccs_file_util.file_exists(shape_dir+"/points.ob") and int(do_refresh)==0):
      sp.ob_file=shape_dir+"/points.ob"
    else:
      agr = Gr()
      [year, day, ts] = model_util.corie2ydn(start_time)
      dname = "%s/%d-%03d/run/" % (base_dir+"/"+estuaryid+"/", year, day)      
      agr.readHGrid(dname+"hgrid.gr3")
      aob = Ob()
      sppts = sp.get_pts("point")
      aob.grid_init(agr, array([sppts[0], sppts[1]]).T)
      aob.save_ob(shape_dir+"/points.ob")
      sp.ob_file = shape_dir+"/points.ob"

    if (mean(abs(points[0][:,2]))==0 and mean(abs(points[0][:,3]))==0) or int(do_refresh)==1:
      sp.gen_channel_norms(base_dir+"/"+estuaryid+"/", start_time, start_time+2, "point", 16, 1, 0)
      meccs_file_util.write_cs_file(shape_dir+"/points.txt", sp.shapes, pnames)      

    mr_pt = mreader(base_dir+"/"+estuaryid+"/", apoints[:,0:2], start_time, end_time, 0, 16, 1, 0, shape_dir+"/points.ob", 1)    

    for i in range(ndays):
      curr_time = start_time+i

      data = sp.gen_point_data(base_dir+"/"+estuaryid+"/", curr_time, curr_time+(1-0.25/24), 16, "base", 1, 0)
      data["start_time"] = curr_time
      data["end_time"] = curr_time+(1-0.25/24)
      data["dt"] = 900      
      dname = "%s/%s/%d-%03d/process/" % (base_dir, estuaryid, yr, day+i)
      sp.save_cs_data(data, dname+"/point_data.dat")

      ##### get raw point data
      (salt, eta, dry, z) = mr_pt.get_model_data("salt.63", curr_time, curr_time+(1-0.25/24))      
      (temp, eta, dry, z) = mr_pt.get_model_data("temp.63", curr_time, curr_time+(1-0.25/24))
      (hvel, eta, dry, z) = mr_pt.get_model_data("hvel.64", curr_time, curr_time+(1-0.25/24))      
      (zcor, eta, dry, z) = mr_pt.get_model_data("zcor.63", curr_time, curr_time+(1-0.25/24))

      for j in range(len(pnames)):
        meccs_file_util.write_raw_data(dname, pnames[j]+"_salt.dat", salt[j,:,0,:], [pnames[j]], curr_time, 900, "seconds", "psu")
        meccs_file_util.write_raw_data(dname, pnames[j]+"_temp.dat", temp[j,:,0,:], [pnames[j]], curr_time, 900, "seconds", "C")
        meccs_file_util.write_raw_data(dname, pnames[j]+"_hvel_u.dat", hvel[j,:,0,:], [pnames[j]], curr_time, 900, "seconds", "m/s")
        meccs_file_util.write_raw_data(dname, pnames[j]+"_hvel_v.dat", hvel[j,:,1,:], [pnames[j]], curr_time, 900, "seconds", "m/s")
        meccs_file_util.write_raw_data(dname, pnames[j]+"_elev.dat", eta[j,:], [pnames[j]], curr_time, 900, "seconds", "m from MSL")
        meccs_file_util.write_raw_data(dname, pnames[j]+"_zmsl.dat", z[j,:], [pnames[j]], curr_time, 900, "seconds", "m from MSL")
        meccs_file_util.write_raw_data(dname, pnames[j]+"_zcor.dat", zcor[j,:,0,:], [pnames[j]], curr_time, 900, "seconds", "m from MSL")
#regions
if (os.path.exists(shape_dir+"/regions.txt")):
  sp = shape_procs.shape_procs()

  (regions, rnames, rrivers) = meccs_file_util.read_cs_file(shape_dir+"/regions.txt", 0)
  
  if len(regions)>0:
    for i in range(len(regions)):
      sp.add_shape(regions[i], "volume", 0, rnames[i])

#    data = sp.gen_volume_data(base_dir+"/"+estuaryid+"/", start_time, end_time, 16, "set_1", 0)
#    data["start_time"] = start_time
#    data["end_time"] = end_time
#    sp.save_cs_data(data, out_base+"/volume_data_"+str(int(start_time))+"_"+str(int(end_time))+".dat")

    for i in range(ndays):
      curr_time = start_time+i

      data = sp.gen_volume_data(base_dir+"/"+estuaryid+"/", curr_time, curr_time+(1-0.25/24), 16, "base", 0)
      data["start_time"] = curr_time
      data["end_time"] = curr_time+(1-0.25/24)
      data["dt"] = 900      
      dname = "%s/%s/%d-%03d/process/" % (base_dir, estuaryid, yr, day+i)
      sp.save_cs_data(data, dname+"/volume_data.dat")
