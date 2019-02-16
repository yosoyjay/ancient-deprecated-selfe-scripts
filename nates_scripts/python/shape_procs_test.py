from numpy import *
from gr import *
import model_util
from model_reader_test import mreader_test
from geometry import *
#import pdb

# shape process class calculates physical properties accross various shapes:
#   supported shapes will eventually be: point, cross section, transect and volume
#   current (1/2/2008) status is that some cross section

class shape_procs_test:
  class _struct:
    pass  
  def __init__(self):
    self.types_enum = array(["point", "cross_sec", "transect", "volume"])
    self.ob_file = 0
 
  def add_shape(self, shape, stype, res, name):
    #shape is an n X 4 array like [x pts, y pts, x normal, y normal]
    #res is the resolution of the line:, e.g. points every 50 m
    #  set to 0 to not resolve the line
    #  resolved normals are the average of their endpoint normals (normal vectors for along-channel projection)
    #    set to 0 if projection is not desired
    
    tmp = stype == self.types_enum

    if not(tmp.any()):
      print "error, shape type must be one of "+str(self.types_enum)
    else:
      try:
        self.shapes_types = append(self.shapes_types, stype)
        if (res and (stype == "cross_sec" or stype == "transect")):
          self.shapes.append(model_util.resolve_line(shape[:,0], shape[:,1], shape[:,2], shape[:,3], res))
        else:
          self.shapes.append(shape)
        self.names = append(self.names, name)
      except AttributeError:
        self.shapes_types = array([stype])
        if (res and (stype == "cross_sec" or stype == "transect")):
          self.shapes = [model_util.resolve_line(shape[:,0], shape[:,1], shape[:,2], shape[:,3], res)]
        else:
          self.shapes = [shape]
        self.names = array([name])

  def get_shape(self, name):
    ret = nan
    indx = where(self.names==name)
    if (len(indx)>0):
      ret = self.shapes[indx[0][0]]
    return ret

  def get_pts(self, shape_type, do_norms=0):
    indx = where(self.shapes_types == shape_type)[0]
    xpts = array([])
    ypts = array([])
    xnorms = array([])
    ynorms = array([])
    pts_indx = array([], int)
    names = []
    for i in indx:
      names.append(self.names[i])
      if pts_indx.shape[0]==0:
        pts_indx = append(pts_indx, 0)
      else:
        pts_indx = append(pts_indx, pts_indx[-1]+curr_shape.shape[0])
      curr_shape = self.shapes[i]

      xpts = append(xpts, curr_shape[:,0])
      ypts = append(ypts, curr_shape[:,1])      
      xnorms = append(xnorms, curr_shape[:,2])
      ynorms = append(ynorms, curr_shape[:,3])        
    pts_indx = append(pts_indx, pts_indx[-1]+curr_shape.shape[0])

    return (xpts, ypts, xnorms, ynorms, names, pts_indx)

  def gen_point_data(self, base_dir, start_time, end_time, db, data_type="base", do_norms=0, hindcast=1):
    ret = {}
    if hindcast==1:
      nddir = 7
    else:
      nddir = 1
      
    indx = where(self.shapes_types == "point")[0]
    (xpts, ypts, xnorms, ynorms, names, nada) = self.get_pts("point")
#    indx = where(self.shapes_types == "point")[0]
#    xpts = array([])
#    ypts = array([])
#    xnorms = array([])
#    ynorms = array([])
#    names = []
#    for i in indx:
#      names.append(self.names[i])
#      curr_shape = self.shapes[i]
#
#      xpts = append(xpts, curr_shape[0,0])
#      ypts = append(ypts, curr_shape[0,1])
#      
#      if (do_norms):
#        xnorms = append(xnorms, curr_shape[0,2])
#        ynorms = append(ynorms, curr_shape[0,3])        
       
    if (data_type == "base"):
      mr = self.get_model_reader(base_dir, xpts, ypts, start_time, end_time, 0, db, nddir, hindcast)

        
      (sdat, eta, dry) = mr.get_model_data(dtype="salt.63")
      depths = dot(mr.ob.H, mr.grid.depth)
      z = mr.computeZlevels(depths, zeros((mr.ob.H.shape[0])))
      (tdat, eta, dry) = mr.get_model_data(dtype="temp.63")
      rho = model_util.density(sdat, tdat)              
      (vdat, eta, dry) = mr.get_model_data(dtype="hvel.64")

      bindx = zeros(eta.shape[0], int)
      for i in range(depths.shape[0]):             #should be a clever way to do this without a loop... not right to assume S-level
        tmp = where((z[i,:]+depths[i])>=0)[0]
        if tmp.shape[0] > 0:          
          bindx[i] = tmp.min()
        else:
          bindx[i] = mr.grid.vgrid.kz-1

      davg_salt = zeros((indx.shape[0], sdat.shape[3]), float)
      surf_salt = zeros((indx.shape[0], sdat.shape[3]), float)
      bott_salt = zeros((indx.shape[0], sdat.shape[3]), float)
      davg_den = zeros((indx.shape[0], sdat.shape[3]), float)
      surf_den = zeros((indx.shape[0], sdat.shape[3]), float)
      bott_den = zeros((indx.shape[0], sdat.shape[3]), float)
      davg_temp = zeros((indx.shape[0], sdat.shape[3]), float)
      surf_temp = zeros((indx.shape[0], sdat.shape[3]), float)
      bott_temp = zeros((indx.shape[0], sdat.shape[3]), float)
      davg_hvel_mag = zeros((indx.shape[0], sdat.shape[3]), float)
      surf_hvel_mag = zeros((indx.shape[0], sdat.shape[3]), float)
      davg_hvel_dir = zeros((indx.shape[0], sdat.shape[3]), float)
      surf_hvel_dir = zeros((indx.shape[0], sdat.shape[3]), float)
      depths_ts = zeros(eta.shape, float)

      for i in range(eta.shape[1]):
        depths_ts[:,i] = depths + eta[:,i]
        z = mr.computeZlevels(depths, eta[:,i])

        #deal with dry nodes:
        dindx = where(dry[:,i] == 1)[0]
        depths_ts[dindx, i] = 0

        #widely used calculations:
        for k in range(indx.shape[0]):
          vSegLens = z[k, bindx[k]+1:]-z[k, bindx[k]:-1]
          dry1 = where(dindx==k)[0]

          #get salinity data: cross section avg, cross section normalized surface and bottom
          if dry1.shape[0]>0:
            davg_salt[k,i] = 0
            surf_salt[k,i] = 0
            bott_salt[k,i] = 0
          else:
            davg_salt[k,i] = sum(((sdat[k, bindx[k]:-1, 0, i] + sdat[k, bindx[k]+1:, 0, i])/2)*vSegLens)/depths_ts[k,i]
            surf_salt[k,i] = sdat[k,-1,0,i]
            bott_salt[k,i] = sdat[k,bindx[k],0,i]
          #density
          if dry1.shape[0]>0:
            davg_den[k,i] = 0
            surf_den[k,i] = 0
            bott_den[k,i] = 0
          else:
            davg_den[k,i] = sum(((rho[k, bindx[k]:-1, 0, i] + rho[k, bindx[k]+1:, 0, i])/2)*vSegLens)/depths_ts[k,i]
            surf_den[k,i] = rho[k,-1,0,i]
            bott_den[k,i] = rho[k,bindx[k],0,i]
          #temp
          if dry1.shape[0]>0:
            davg_temp[k,i] = 0
            surf_temp[k,i] = 0
            bott_temp[k,i] = 0
          else:
            davg_temp[k,i] = sum(((tdat[k, bindx[k]:-1, 0, i] + tdat[k, bindx[k]+1:, 0, i])/2)*vSegLens)/depths_ts[k,i]
            surf_temp[k,i] = tdat[k,-1,0,i]
            bott_temp[k,i] = tdat[k,bindx[k],0,i]

          #get velocity data, including magnitude and along-channel.
          if dry1.shape[0]>0:
            davg_hvel_mag[k,i] = 0
            surf_hvel_mag[k,i] = 0
          else:
            davg_hvel_mag[k,i] = sum((magnitude(vdat[k, bindx[k]:-1, 0, i], vdat[k, bindx[k]:-1, 1, i]) + \
                                          magnitude(vdat[k, bindx[k]+1:, 0, i], vdat[k, bindx[k]+1:, 1, i]))/2*vSegLens)/depths_ts[k,i]
            surf_hvel_mag[k,i] = magnitude(vdat[k,-1,0,i], vdat[k,-1,1,i])

          if dry1.shape[0]>0:
            davg_hvel_dir[k,i] = 0
            surf_hvel_dir[k,i] = 0
          else:
            davg_hvel_dir[k,i] = sum((vdat[k, bindx[k]:-1, 0, i]*xnorms[k] + vdat[k, bindx[k]:-1, 1, i]*ynorms[k] + \
                                          vdat[k, bindx[k]+1:, 0, i]*xnorms[k] + vdat[k, bindx[k]+1:, 1, i]*ynorms[k])/2*vSegLens)/depths_ts[k,i]
            surf_hvel_dir[k,i] = vdat[k, -1, 0, i]*xnorms[k] + vdat[k, -1, 1, i]*ynorms[k]

      ret["elev"] = eta
      ret["surf_salt"] = surf_salt
      ret["bott_salt"] = bott_salt
      ret["davg_salt"] = davg_salt      
      ret["davg_den"] = davg_den
      ret["surf_den"] = surf_den
      ret["bott_den"] = bott_den
      ret["davg_temp"] = davg_temp
      ret["surf_temp"] = surf_temp
      ret["bott_temp"] = bott_temp
      ret["davg_hvel_mag"] = davg_hvel_mag
      ret["surf_hvel_mag"] = surf_hvel_mag
      ret["davg_hvel_dir"] = davg_hvel_dir
      ret["surf_hvel_dir"] = surf_hvel_dir
      ret["depths"] = depths
      ret["names"] = names
      ret["dt"] = mr.dt
    
    return ret

  def gen_transect_data(self, base_dir, start_time, end_time, name, db, do_norms=0, ndays_dir=7, hindcast=1):
    ret = {}

    tmp = where(self.names == name)

    if (len(tmp)>0):
      indx = tmp[0][0]
      names = self.names[indx]
      curr_shape = self.shapes[indx]
  
      xpts = curr_shape[:,0]
      ypts = curr_shape[:,1]
      
#      if (do_norms):
      xnorms = curr_shape[:,2]
      ynorms = curr_shape[:,3]

      mr = self.get_model_reader(base_dir, xpts, ypts, start_time, end_time, 0, db, ndays_dir, hindcast)

      (sdat, eta, dry) = mr.get_model_data(dtype="salt.63")
      depths = dot(mr.ob.H, mr.grid.depth)
      z = mr.computeZlevels(depths, zeros((mr.ob.H.shape[0])))
      (tdat, eta, dry) = mr.get_model_data(dtype="temp.63")
      rho = model_util.density(sdat, tdat)
      (vdat, eta, dry) = mr.get_model_data(dtype="hvel.64")
  
      bindx = zeros(eta.shape[0], int)
      for i in range(depths.shape[0]):             #should be a clever way to do this without a loop... not right to assume S-level
        tmp = where((z[i,:]+depths[i])>=0)[0]
        if tmp.shape[0] > 0:          
          bindx[i] = tmp.min()
        else:
          bindx[i] = mr.grid.vgrid.kz-1

      dpth = zeros((sdat.shape[0], sdat.shape[3]), float)
      dmin_salt = zeros((sdat.shape[0], sdat.shape[3]), float)
      dmax_salt = zeros((sdat.shape[0], sdat.shape[3]), float)
      surf_salt = zeros((sdat.shape[0], sdat.shape[3]), float)
      bott_salt = zeros((sdat.shape[0], sdat.shape[3]), float)
      davg_salt = zeros((sdat.shape[0], sdat.shape[3]), float)
      dmin_den = zeros((sdat.shape[0], sdat.shape[3]), float)
      dmax_den = zeros((sdat.shape[0], sdat.shape[3]), float)    
      davg_den = zeros((sdat.shape[0], sdat.shape[3]), float)
      surf_den = zeros((sdat.shape[0], sdat.shape[3]), float)
      bott_den = zeros((sdat.shape[0], sdat.shape[3]), float)
      dmin_temp = zeros((sdat.shape[0], sdat.shape[3]), float)
      dmax_temp = zeros((sdat.shape[0], sdat.shape[3]), float)
      davg_temp = zeros((sdat.shape[0], sdat.shape[3]), float)
      surf_temp = zeros((sdat.shape[0], sdat.shape[3]), float)
      bott_temp = zeros((sdat.shape[0], sdat.shape[3]), float)
      dmin_hvel_mag = zeros((sdat.shape[0], sdat.shape[3]), float)
      dmax_hvel_mag = zeros((sdat.shape[0], sdat.shape[3]), float)    
      davg_hvel_mag = zeros((sdat.shape[0], sdat.shape[3]), float)
      surf_hvel_mag = zeros((sdat.shape[0], sdat.shape[3]), float)
      dmin_hvel_dir = zeros((sdat.shape[0], sdat.shape[3]), float)
      dmax_hvel_dir = zeros((sdat.shape[0], sdat.shape[3]), float)    
      davg_hvel_dir = zeros((sdat.shape[0], sdat.shape[3]), float)
      surf_hvel_dir = zeros((sdat.shape[0], sdat.shape[3]), float)    
      depths_ts = zeros(eta.shape, float)

      for i in range(eta.shape[1]):        
        depths_ts[:,i] = depths + eta[:,i]
        z = mr.computeZlevels(depths, eta[:,i])
  
        #deal with dry nodes:
        idry = where(dry[:,i] == 1)[0]
        depths_ts[idry, i] = 0

#        if i%10 == 0:
#          pdb.set_trace()
  
        for j in range(xpts.shape[0]):
  
          #widely used calculations:
          vSegLens1 = z[j, bindx[j]+1:]-z[j, bindx[j]:-1]
          dry1 = where(idry==j)[0]
  
          if (dry1.shape[0]>0):
            dpth[j,i]          = nan
            dmin_salt[j,i]     = nan
            dmax_salt[j,i]     = nan
            davg_salt[j,i]     = nan
            surf_salt[j,i]     = nan
            bott_salt[j,i]     = nan
            dmin_den[j,i]      = nan
            dmax_den[j,i]      = nan
            davg_den[j,i]      = nan
            surf_den[j,i]      = nan
            bott_den[j,i]      = nan
            dmin_temp[j,i]     = nan
            dmax_temp[j,i]     = nan
            davg_temp[j,i]     = nan
            surf_temp[j,i]     = nan
            bott_temp[j,i]     = nan
            dmin_hvel_mag[j,i] = nan
            dmax_hvel_mag[j,i] = nan
            dmin_hvel_dir[j,i] = nan
            dmax_hvel_dir[j,i] = nan
            davg_hvel_mag[j,i] = nan
            surf_hvel_mag[j,i] = nan
            davg_hvel_dir[j,i] = nan
            surf_hvel_dir[j,i] = nan
          else:
            dpth[j,i]          = depths_ts[j,i]
            dmin_salt[j,i]     = min(sdat[j, bindx[j]:, 0, i])
            dmax_salt[j,i]     = max(sdat[j, bindx[j]:, 0, i])            
            davg_salt[j,i]     = sum(((sdat[j, bindx[j]:-1, 0, i] + sdat[j, bindx[j]+1:, 0, i])/2)*vSegLens1)/sum(vSegLens1)
            surf_salt[j,i]     = sdat[j,-1,0,i]      
            bott_salt[j,i]     = sdat[j,bindx[j],0,i]
            dmin_den[j,i]      = min(rho[j, bindx[j]:, 0, i])
            dmax_den[j,i]      = max(rho[j, bindx[j]:, 0, i])
            davg_den[j,i]      = sum(((rho[j, bindx[j]:-1, 0, i] + rho[j, bindx[j]+1:, 0, i])/2)*vSegLens1)/sum(vSegLens1)
            surf_den[j,i]      = rho[j,-1,0,i]
            bott_den[j,i]      = rho[j,bindx[j],0,i]
            dmin_temp[j,i]     = min(tdat[j, bindx[j]:, 0, i])
            dmax_temp[j,i]     = max(tdat[j, bindx[j]:, 0, i])
            davg_temp[j,i]     = sum(((tdat[j, bindx[j]:-1, 0, i] + tdat[j, bindx[j]+1:, 0, i])/2)*vSegLens1)/sum(vSegLens1)
            surf_temp[j,i]     = tdat[j,-1,0,i]      
            bott_temp[j,i]     = tdat[j,bindx[j],0,i]
            dmin_hvel_mag[j,i] = min(magnitude(vdat[j, bindx[j]+1:, 0, i], vdat[j,bindx[j]+1:, 1, i]))            
            dmax_hvel_mag[j,i] = max(magnitude(vdat[j, bindx[j]:, 0, i], vdat[j,bindx[j]:, 1, i]))
            dmin_hvel_dir[j,i] = min(vdat[j, bindx[j]+1:, 0, i]*xnorms[j] + vdat[j,bindx[j]+1:, 1, i]*ynorms[j])
            dmax_hvel_dir[j,i] = max(vdat[j, bindx[j]:, 0, i]*xnorms[j] + vdat[j,bindx[j]:, 1, i]*ynorms[j])
            davg_hvel_mag[j,i] = sum((magnitude(vdat[j, bindx[j]:-1, 0, i], vdat[j, bindx[j]:-1, 1, i]) + \
                      magnitude(vdat[j, bindx[j]+1:, 0, i], vdat[j, bindx[j]+1:, 1, i]))/2*vSegLens1)/sum(vSegLens1)
            surf_hvel_mag[j,i] = magnitude(vdat[j,-1,0,i], vdat[j,-1,1,i])
            davg_hvel_dir[j,i] = sum(((vdat[j, bindx[j]:-1, 0, i]*xnorms[j] + vdat[j, bindx[j]:-1, 1, i]*ynorms[j]) + \
                        (vdat[j, bindx[j]+1:, 0, i]*xnorms[j] + vdat[j, bindx[j]+1:, 1, i]*ynorms[j]))/2*vSegLens1)/sum(vSegLens1)
            surf_hvel_dir[j,i] = vdat[j,-1,0,i]*xnorms[j]+vdat[j,-1,1,i]*ynorms[j]


#      pdb.set_trace()
            
      ret["depth"] = dpth
      ret["dmin_salt"] = dmin_salt 
      ret["dmax_salt"] = dmax_salt
      ret["surf_salt"] = surf_salt
      ret["bott_salt"] = bott_salt
      ret["davg_salt"] = davg_salt
      ret["dmin_den"] = dmin_den
      ret["dmax_den"] = dmax_den
      ret["davg_den"] = davg_den
      ret["surf_den"] = surf_den
      ret["bott_den"] = bott_den
      ret["dmin_temp"] = dmin_temp
      ret["dmax_temp"] = dmax_temp
      ret["davg_temp"] = davg_temp
      ret["surf_temp"] = surf_temp
      ret["bott_temp"] = bott_temp
      ret["dmin_hvel_mag"] = dmin_hvel_mag
      ret["dmax_hvel_mag"] = dmax_hvel_mag 
      ret["davg_hvel_mag"] = davg_hvel_mag
      ret["surf_hvel_mag"] = surf_hvel_mag
      ret["dmin_hvel_dir"] = dmin_hvel_dir
      ret["dmax_hvel_dir"] = dmax_hvel_dir 
      ret["davg_hvel_dir"] = davg_hvel_dir
      ret["surf_hvel_dir"] = surf_hvel_dir
      ret["dt"] = mr.dt
      ret["name"] = name
        
    return ret    

  def gen_cross_section_data(self, base_dir, start_time, end_time, db, data_type="base", do_norms=0, ndays_dir=7, hindcast=1):
    #collect the points in the cross sections:
    ret = {}

    (xpts, ypts, xnorms, ynorms, names, pts_indx) = self.get_pts("cross_sec")
    indx = where(self.shapes_types == "cross_sec")[0]

#    xpts = array([])
#    ypts = array([])
#    xnorms = array([])
#    ynorms = array([])
#    pts_indx = array([], int)
#    names = []
#    for i in indx:
#      names.append(self.names[i])
#      if pts_indx.shape[0]==0:
#        pts_indx = append(pts_indx, 0)
#      else:
#        pts_indx = append(pts_indx, pts_indx[-1]+curr_shape.shape[0])
#      curr_shape = self.shapes[i]
#
#      xpts = append(xpts, curr_shape[:,0])
#      ypts = append(ypts, curr_shape[:,1])
#      
#      if (do_norms):
#        xnorms = append(xnorms, curr_shape[:,2])
#        ynorms = append(ynorms, curr_shape[:,3])        
#    pts_indx = append(pts_indx, pts_indx[-1]+curr_shape.shape[0])

    if data_type == "base":     #a basic package of: cs area, depth/surface/bottom avg salt,
                                #hvel (along-channel and magnitude), temp, density

      mr = self.get_model_reader(base_dir, xpts, ypts, start_time, end_time, 0, db, ndays_dir, hindcast)
        
      (sdat, eta, dry) = mr.get_model_data(dtype="salt.63")
      depths = dot(mr.ob.H, mr.grid.depth)
      z = mr.computeZlevels(depths, zeros((mr.ob.H.shape[0])))
      (tdat, eta, dry) = mr.get_model_data(dtype="temp.63")
      rho = model_util.density(sdat, tdat)
      (vdat, eta, dry) = mr.get_model_data(dtype="hvel.64")

      bindx = zeros(eta.shape[0], int)
      for i in range(depths.shape[0]):             #should be a clever way to do this without a loop... not right to assume S-level
        tmp = where((z[i,:]+depths[i])>=0)[0]
        if tmp.shape[0] > 0:          
          bindx[i] = tmp.min()
        else:
          bindx[i] = mr.grid.vgrid.kz-1

      cs_areas = zeros((indx.shape[0], sdat.shape[3]), float)
      davg_salt = zeros((indx.shape[0], sdat.shape[3]), float)
      cs_surf_salt = zeros((indx.shape[0], sdat.shape[3]), float)
      cs_bott_salt = zeros((indx.shape[0], sdat.shape[3]), float)
      davg_den = zeros((indx.shape[0], sdat.shape[3]), float)
      cs_surf_den = zeros((indx.shape[0], sdat.shape[3]), float)
      cs_bott_den = zeros((indx.shape[0], sdat.shape[3]), float)
      davg_temp = zeros((indx.shape[0], sdat.shape[3]), float)
      cs_surf_temp = zeros((indx.shape[0], sdat.shape[3]), float)
      cs_bott_temp = zeros((indx.shape[0], sdat.shape[3]), float)
      davg_hvel_mag = zeros((indx.shape[0], sdat.shape[3]), float)
      cs_surf_hvel_mag = zeros((indx.shape[0], sdat.shape[3]), float)
      davg_hvel_dir = zeros((indx.shape[0], sdat.shape[3]), float)
      cs_surf_hvel_dir = zeros((indx.shape[0], sdat.shape[3]), float)
      hvel_flux = zeros((indx.shape[0], sdat.shape[3]), float)
      cs_depth = zeros((indx.shape[0], sdat.shape[3]), float)
      
      depths_ts = zeros(eta.shape, float)
      for i in range(eta.shape[1]):
        depths_ts[:,i] = depths + eta[:,i]
        z = mr.computeZlevels(depths, eta[:,i])

        #deal with dry nodes:
        idry = where(dry[:,i] == 1)[0]
        depths_ts[idry, i] = 0

        #get the cross sectional areas
        for j in range(pts_indx.shape[0]-1):
          for k in range(pts_indx[j], pts_indx[j+1]-1):
            tSegLen = distance(xpts[k], ypts[k], xpts[k+1], ypts[k+1])

            #widely used calculations:
            vSegLens1 = z[k, bindx[k]+1:]-z[k, bindx[k]:-1]
            vSegLens2 = z[k+1, bindx[k+1]+1:]-z[k+1, bindx[k+1]:-1]
            dry1 = where(idry==k)[0]
            dry2 = where(idry==k+1)[0]            

            #make sure cs area doesn't include dry pts
            hSegLen = tSegLen

            if (dry1.shape[0]>0):
              hSegLen = hSegLen-tSegLen/2.
              depth1 = 0
              salt1 = 0
              den1 = 0
              sden1 = 0
              bden1 = 0
              temp1 = 0
              hvel1m = 0
              hvel1d = 0              
            else:
              depth1 = depths_ts[k,i]                
              salt1 = sum(((sdat[k, bindx[k]:-1, 0, i] + sdat[k, bindx[k]+1:, 0, i])/2)*vSegLens1)
              den1 = sum(((rho[k, bindx[k]:-1, 0, i] + rho[k, bindx[k]+1:, 0, i])/2)*vSegLens1)
              sden1 = rho[k,-1,0,i]*depths_ts[k,i]
              bden1 = rho[k,bindx[k],0,i]*depths_ts[k,i]
              temp1 = sum(((tdat[k, bindx[k]:-1, 0, i] + tdat[k, bindx[k]+1:, 0, i])/2)*vSegLens1)
              hvel1m = (magnitude(vdat[k, bindx[k]:-1, 0, i], vdat[k, bindx[k]:-1, 1, i]) + \
                        magnitude(vdat[k, bindx[k]+1:, 0, i], vdat[k, bindx[k]+1:, 1, i]))/2*vSegLens1              
              hvel1d = ((vdat[k, bindx[k]:-1, 0, i]*xnorms[k] + vdat[k, bindx[k]:-1, 1, i]*ynorms[k]) + \
                        (vdat[k, bindx[k]+1:, 0, i]*xnorms[k] + vdat[k, bindx[k]+1:, 1, i]*ynorms[k]))/2*vSegLens1              
            if (dry2.shape[0]>0):
              hSegLen = hSegLen-tSegLen/2.
              depth2 = 0
              salt2 = 0
              den2 = 0
              sden2 = 0
              bden2 = 0 
              temp2 = 0
              hvel2m = 0
              hvel2d = 0              
            else:
              depth2 = depths_ts[k+1,i]
              salt2 = sum(((sdat[k+1, bindx[k+1]:-1, 0, i] + sdat[k+1, bindx[k+1]+1:, 0, i])/2)*vSegLens2)
              den2 = sum(((rho[k+1, bindx[k+1]:-1, 0, i] + rho[k+1, bindx[k+1]+1:, 0, i])/2)*vSegLens2)
              sden2 = rho[k+1,-1,0,i]*depths_ts[k+1,i]
              bden2 = rho[k+1,bindx[k+1],0,i]*depths_ts[k+1,i]
              temp2 = sum(((tdat[k+1, bindx[k+1]:-1, 0, i] + tdat[k+1, bindx[k+1]+1:, 0, i])/2)*vSegLens2)
              hvel2m = (magnitude(vdat[k+1, bindx[k+1]:-1, 0, i], vdat[k+1, bindx[k+1]:-1, 1, i]) + \
                        magnitude(vdat[k+1, bindx[k+1]+1:, 0, i], vdat[k+1, bindx[k+1]+1:, 1, i]))/2*vSegLens2
              hvel2d = ((vdat[k+1, bindx[k+1]:-1, 0, i]*xnorms[k+1] + vdat[k+1, bindx[k+1]:-1, 1, i]*ynorms[k]) + \
                       (vdat[k+1, bindx[k+1]+1:, 0, i]*xnorms[k+1] + vdat[k+1, bindx[k+1]+1:, 1, i]*ynorms[k]))/2*vSegLens2              
            cs_areas[j,i] = cs_areas[j,i] + (depth1 + depth2)/2 * hSegLen
            cs_depth[j,i] = cs_depth[j,i] + (depth1**2 + depth2**2)/2 * hSegLen

            #get salinity data: cross section avg, cross section normalized surface and bottom

#            davg_salt[j,i] = davg_salt[j,i]+sum((salt1+salt2)*hSegLen/2)
            davg_salt[j,i] = davg_salt[j,i]+(salt1+salt2)*hSegLen/2
            cs_surf_salt[j,i] = cs_surf_salt[j,i]+(sdat[k,-1,0,i]*depths_ts[k,i] + sdat[k+1,-1,0,i]*depths_ts[k+1,i])*hSegLen/2
            cs_bott_salt[j,i] = cs_bott_salt[j,i]+(sdat[k,bindx[k],0,i]*depths_ts[k,i] + sdat[k+1,bindx[k+1],0,i]*depths_ts[k+1,i])*hSegLen/2

            #get density data: (same as salt)
            #density() returns nan when salt or temp are negative (ie are dry), so surf and bott den need to account for nan

            davg_den[j,i] = davg_den[j,i]+(den1+den2)*hSegLen/2
            cs_surf_den[j,i] = cs_surf_den[j,i]+(sden1 + sden2)*hSegLen/2
            cs_bott_den[j,i] = cs_bott_den[j,i]+(bden1 + bden2)*hSegLen/2

            #get temp data:
            davg_temp[j,i] = davg_temp[j,i]+(temp1+temp2)*hSegLen
            cs_surf_temp[j,i] = cs_surf_temp[j,i]+(tdat[k,-1,0,i]*depths_ts[k,i] + tdat[k+1,-1,0,i]*depths_ts[k+1,i])*hSegLen/2
            cs_bott_temp[j,i] = cs_bott_temp[j,i]+(tdat[k,bindx[k],0,i]*depths_ts[k,i] + tdat[k+1,bindx[k+1],0,i]*depths_ts[k+1,i])*hSegLen/2

            #get velocity data, including magnitude and along-channel.  Calculate cross sectional flux as well:
 
            davg_hvel_mag[j,i] = davg_hvel_mag[j,i]+sum((hvel1m+hvel2m)*hSegLen/2)
            cs_surf_hvel_mag[j,i] = cs_surf_hvel_mag[j,i]+(magnitude(vdat[k,-1,0,i], vdat[k,-1,1,i])*depths_ts[k,i] + \
                                                           magnitude(vdat[k+1,-1,0,i], vdat[k+1,-1,1,i])*depths_ts[k+1,i])*hSegLen/2
            davg_hvel_dir[j,i] = davg_hvel_dir[j,i]+sum((hvel1d+hvel2d)*hSegLen/2)
            cs_surf_hvel_dir[j,i] = cs_surf_hvel_dir[j,i]+((vdat[k,-1,0,i]*xnorms[k]+vdat[k,-1,1,i]*ynorms[k])*depths_ts[k,i] + \
                                                           (vdat[k+1,-1,0,i]*xnorms[k+1]+vdat[k+1,-1,1,i]*ynorms[k+1])*depths_ts[k+1,i])*hSegLen/2


        #divide by normalizing area
        davg_salt[:,i] = davg_salt[:,i]/cs_areas[:,i]
        cs_surf_salt[:,i] = cs_surf_salt[:,i]/cs_areas[:,i]
        cs_bott_salt[:,i] = cs_bott_salt[:,i]/cs_areas[:,i]
        cs_depth[:,i] = cs_depth[:,i]/cs_areas[:,i]

        davg_den[:,i] = davg_den[:,i]/cs_areas[:,i]
        cs_surf_den[:,i] = cs_surf_den[:,i]/cs_areas[:,i]
        cs_bott_den[:,i] = cs_bott_den[:,i]/cs_areas[:,i]
        davg_temp[:,i] = davg_temp[:,i]/cs_areas[:,i]
        cs_surf_temp[:,i] = cs_surf_temp[:,i]/cs_areas[:,i]
        cs_bott_temp[:,i] = cs_bott_temp[:,i]/cs_areas[:,i]               
        davg_hvel_mag[:,i] = davg_hvel_mag[:,i]/cs_areas[:,i]
        cs_surf_hvel_mag[:,i] = cs_surf_hvel_mag[:,i]/cs_areas[:,i]
        hvel_flux[:,i] = davg_hvel_dir[:,i]   
        davg_hvel_dir[:,i] = davg_hvel_dir[:,i]/cs_areas[:,i]        
        cs_surf_hvel_dir[:,i] = cs_surf_hvel_dir[:,i]/cs_areas[:,i]

      ret["davg_salt"] = davg_salt
      ret["cs_surf_salt"] = cs_surf_salt
      ret["cs_bott_salt"] = cs_bott_salt
      ret["davg_den"] = davg_den
      ret["cs_surf_den"] = cs_surf_den
      ret["cs_bott_den"] = cs_bott_den
      ret["davg_temp"] = davg_temp
      ret["cs_surf_temp"] = cs_surf_temp
      ret["cs_bott_temp"] = cs_bott_temp
      ret["davg_hvel_mag"] = davg_hvel_mag
      ret["cs_surf_hvel_mag"] = cs_surf_hvel_mag
      ret["hvel_flux"] = hvel_flux
      ret["davg_hvel_dir"] = davg_hvel_dir
      ret["cs_surf_hvel_dir"] = cs_surf_hvel_dir
      ret["names"] = names      
      ret["dt"] = mr.dt
      ret["cs_areas"] = cs_areas
      ret["cs_depth"] = cs_depth      
      
    return ret

  def get_model_reader(self, base_dir, xpts, ypts, start_time, end_time, dvg, db, ndd, hindcast, force_new_ob=0):
    if (self.ob_file and force_new_ob==0):
      mr = mreader_test(base_dir, array([xpts, ypts]).T, start_time, end_time, dvg, ndd, db, hindcast, self.ob_file, 1)
    else:
      mr = mreader_test(bd = base_dir, bp = array([xpts, ypts]).T, start=start_time, end=end_time, dvg=dvg, db=db, ndays_dir=ndd, hindcast=hindcast)
    return mr

  def gen_cross_section_data2(self, base_dir, start_time, end_time, db, data_type, vect=0, ndays_dir=7, hindcast=1, mr=0):
    ret = {}

    #collect the points in the cross sections:

    (xpts, ypts, xnorms, ynorms, names, pts_indx) = self.get_pts("cross_sec")
    
#    indx = where(self.shapes_types == "cross_sec")[0]
#    xpts = array([])
#    ypts = array([])
#    xnorms = array([])
#    ynorms = array([])
#    pts_indx = array([], int)
#    names = []
#    for i in indx:
#      names.append(self.names[i])
#      if pts_indx.shape[0]==0:
#        pts_indx = append(pts_indx, 0)
#      else:
#        pts_indx = append(pts_indx, pts_indx[-1]+curr_shape.shape[0])
#      curr_shape = self.shapes[i]
#
#      xpts = append(xpts, curr_shape[:,0])
#      ypts = append(ypts, curr_shape[:,1])
#
#      if (vect):
#        xnorms = append(xnorms, curr_shape[:,2])
#        ynorms = append(ynorms, curr_shape[:,3])        
#    pts_indx = append(pts_indx, pts_indx[-1]+curr_shape.shape[0])

    if (mr==0):
      mr = self.get_model_reader(base_dir, xpts, ypts, start_time, end_time, 0, db, ndays_dir, hindcast)
#      if (self.ob_file):
#        mr = mreader_test(bd = base_dir, bp = array([xpts, ypts]).T, start=start_time, end=end_time, dvg=0, db=db, hindcast=hindcast, self.ob_file, 1)
#      else:
#        mr = mreader_test(bd = base_dir, bp = array([xpts, ypts]).T, start=start_time, end=end_time, dvg=0, db=db, hindcast=hindcast)
      
    if (data_type=="rho"):
      (sdat, eta, dry) = mr.get_model_data(dtype="salt.63")      
      (tdat, eta, dry) = mr.get_model_data(dtype="temp.63")
      dat = model_util.density(sdat, tdat)
    else:      
      (dat, eta, dry) = mr.get_model_data(dtype=data_type)
    depths = dot(mr.ob.H, mr.grid.depth)
    z = mr.computeZlevels(depths, zeros((mr.ob.H.shape[0])))
    bindx = zeros(eta.shape[0], int)
    for i in range(depths.shape[0]):             #should be a clever way to do this without a loop... not right to assume S-level
      tmp = where((z[i,:]+depths[i])>=0)[0]
      if tmp.shape[0] > 0:          
        bindx[i] = tmp.min()
      else:
        bindx[i] = mr.grid.vgrid.kz-1

    cs_areas = zeros((indx.shape[0], sdat.shape[3]), float)
    davg_dat = zeros((indx.shape[0], sdat.shape[3]), float)
    cs_surf_dat = zeros((indx.shape[0], sdat.shape[3]), float)
    cs_bott_dat = zeros((indx.shape[0], sdat.shape[3]), float)
    if vect:
      davg_dat_dir = zeros((indx.shape[0], sdat.shape[3]), float)
      cs_surf_dat_dir = zeros((indx.shape[0], sdat.shape[3]), float)
      cs_bott_dat_dir = zeros((indx.shape[0], sdat.shape[3]), float)
      flux = zeros((indx.shape[0], sdat.shape[3]), float)        
      
    depths_ts = zeros(eta.shape, float)
    for i in range(eta.shape[1]):
      depths_ts[:,i] = depths + eta[:,i]
      z = mr.computeZlevels(depths, eta[:,i])

      #deal with dry nodes:
      idry = where(dry[:,i] == 1)[0]
      depths_ts[idry, i] = 0

      #get the cross sectional areas
      for j in range(pts_indx.shape[0]-1):
        for k in range(pts_indx[j], pts_indx[j+1]-1):
          tSegLen = distance(xpts[k], ypts[k], xpts[k+1], ypts[k+1])

            #widely used calculations:
          vSegLens1 = z[k, bindx[k]+1:]-z[k, bindx[k]:-1]
          vSegLens2 = z[k+1, bindx[k+1]+1:]-z[k+1, bindx[k+1]:-1]
          dry1 = where(idry==k)[0]
          dry2 = where(idry==k+1)[0]            

          #here: making a function that is more general...
          #make sure cs area doesn't include dry pts
          hSegLen = tSegLen

#          if (dry1.shape[0]>0):
#            hSegLen = hSegLen-tSegLen/2.
#            dat1 = 0
#            if vect
#              hvel1m = 0
#              hvel1d = 0              
#            else:
#              depth1 = depths_ts[k,i]                
#              salt1 = sum(((sdat[k, bindx[k]:-1, 0, i] + sdat[k, bindx[k]+1:, 0, i])/2)*vSegLens1)
#              den1 = sum(((rho[k, bindx[k]:-1, 0, i] + rho[k, bindx[k]+1:, 0, i])/2)*vSegLens1)
#              sden1 = rho[k,-1,0,i]*depths_ts[k,i]
#              bden1 = rho[k,bindx[k],0,i]*depths_ts[k,i]
#              temp1 = sum((tdat[k, bindx[k]:-1, 0, i] + tdat[k, bindx[k]+1:, 0, i])/2)*vSegLens1)
#              hvel1m = (magnitude(vdat[k, bindx[k]:-1, 0, i], vdat[k, bindx[k]:-1, 1, i]) + \
#                        magnitude(vdat[k, bindx[k]+1:, 0, i], vdat[k, bindx[k]+1:, 1, i]))/2*vSegLens1              
#              hvel1d = ((vdat[k, bindx[k]:-1, 0, i]*xnorms[k] + vdat[k, bindx[k]:-1, 1, i]*ynorms[k]) + \
#                        (vdat[k, bindx[k]+1:, 0, i]*xnorms[k] + vdat[k, bindx[k]+1:, 1, i]*ynorms[k]))/2*vSegLens1              
#            if (dry2.shape[0]>0):
#              hSegLen = hSegLen-tSegLen/2.
#              depth2 = 0
#              salt2 = 0
#              den2 = 0
#              sden2 = 0
#              bden2 = 0 
#              temp2 = 0
#              hvel2m = 0
#              hvel2d = 0              
#            else:
#              depth2 = depths_ts[k+1,i]
#              salt2 = sum(((sdat[k+1, bindx[k+1]:-1, 0, i] + sdat[k+1, bindx[k+1]+1:, 0, i])/2)*vSegLens2)
#              den2 = sum(((rho[k+1, bindx[k+1]:-1, 0, i] + rho[k+1, bindx[k+1]+1:, 0, i])/2)*vSegLens2)
#              sden2 = rho[k+1,-1,0,i]*depths_ts[k+1,i]
#              bden2 = rho[k+1,bindx[k+1],0,i]*depths_ts[k+1,i]
#              temp2 = sum(((tdat[k+1, bindx[k+1]:-1, 0, i] + tdat[k+1, bindx[k+1]+1:, 0, i])/2)*vSegLens2)
#              hvel2m = (magnitude(vdat[k+1, bindx[k+1]:-1, 0, i], vdat[k+1, bindx[k+1]:-1, 1, i]) + \
#                        magnitude(vdat[k+1, bindx[k+1]+1:, 0, i], vdat[k+1, bindx[k+1]+1:, 1, i]))/2*vSegLens2
#              hvel2d = ((vdat[k+1, bindx[k+1]:-1, 0, i]*xnorms[k+1] + vdat[k+1, bindx[k+1]:-1, 1, i]*ynorms[k]) + \
#                       (vdat[k+1, bindx[k+1]+1:, 0, i]*xnorms[k+1] + vdat[k+1, bindx[k+1]+1:, 1, i]*ynorms[k]))/2*vSegLens2              
#            cs_areas[j,i] = cs_areas[j,i] + (depth1 + depth2)/2 * hSegLen
#
#            #get salinity data: cross section avg, cross section normalized surface and bottom
#
##            davg_salt[j,i] = davg_salt[j,i]+sum((salt1+salt2)*hSegLen/2)
#            davg_salt[j,i] = davg_salt[j,i]+(salt1+salt2)*hSegLen/2
#            cs_surf_salt[j,i] = cs_surf_salt[j,i]+(sdat[k,-1,0,i]*depths_ts[k,i] + sdat[k+1,-1,0,i]*depths_ts[k+1,i])*hSegLen/2
#            cs_bott_salt[j,i] = cs_bott_salt[j,i]+(sdat[k,bindx[k],0,i]*depths_ts[k,i] + sdat[k+1,bindx[k+1],0,i]*depths_ts[k+1,i])*hSegLen/2
#
#            #get density data: (same as salt)
#            #density() returns nan when salt or temp are negative (ie are dry), so surf and bott den need to account for nan
#
#            davg_den[j,i] = davg_den[j,i]+(den1+den2)*hSegLen/2
#            cs_surf_den[j,i] = cs_surf_den[j,i]+(sden1 + sden2)*hSegLen/2
#            cs_bott_den[j,i] = cs_bott_den[j,i]+(bden1 + bden2)*hSegLen/2
#
#            #get temp data:
#            davg_temp[j,i] = davg_temp[j,i]+(temp1+temp2)*hSegLen
#            cs_surf_temp[j,i] = cs_surf_temp[j,i]+(tdat[k,-1,0,i]*depths_ts[k,i] + tdat[k+1,-1,0,i]*depths_ts[k+1,i])*hSegLen/2
#            cs_bott_temp[j,i] = cs_bott_temp[j,i]+(tdat[k,bindx[k],0,i]*depths_ts[k,i] + tdat[k+1,bindx[k+1],0,i]*depths_ts[k+1,i])*hSegLen/2
#
#            #get velocity data, including magnitude and along-channel.  Calculate cross sectional flux as well:
# 
#            davg_hvel_mag[j,i] = davg_hvel_mag[j,i]+sum((hvel1m+hvel2m)*hSegLen/2)
#            cs_surf_hvel_mag[j,i] = cs_surf_hvel_mag[j,i]+(magnitude(vdat[k,-1,0,i], vdat[k,-1,1,i])*depths_ts[k,i] + \
#                                                           magnitude(vdat[k+1,-1,0,i], vdat[k+1,-1,1,i])*depths_ts[k+1,i])*hSegLen/2
#            davg_hvel_dir[j,i] = davg_hvel_dir[j,i]+sum((hvel1d+hvel2d)*hSegLen/2)
#            cs_surf_hvel_dir[j,i] = cs_surf_hvel_dir[j,i]+((vdat[k,-1,0,i]*xnorms[k]+vdat[k,-1,1,i]*ynorms[k])*depths_ts[k,i] + \
#                                                           (vdat[k+1,-1,0,i]*xnorms[k+1]+vdat[k+1,-1,1,i]*ynorms[k+1])*depths_ts[k+1,i])*hSegLen/2
#
#
#        #divide by normalizing area
#        davg_salt[:,i] = davg_salt[:,i]/cs_areas[:,i]
#        cs_surf_salt[:,i] = cs_surf_salt[:,i]/cs_areas[:,i]
#        cs_bott_salt[:,i] = cs_bott_salt[:,i]/cs_areas[:,i]
#
#        davg_den[:,i] = davg_den[:,i]/cs_areas[:,i]
#        cs_surf_den[:,i] = cs_surf_den[:,i]/cs_areas[:,i]
#        cs_bott_den[:,i] = cs_bott_den[:,i]/cs_areas[:,i]
#        davg_temp[:,i] = davg_temp[:,i]/cs_areas[:,i]
#        cs_surf_temp[:,i] = cs_surf_temp[:,i]/cs_areas[:,i]
#        cs_bott_temp[:,i] = cs_bott_temp[:,i]/cs_areas[:,i]               
#        davg_hvel_mag[:,i] = davg_hvel_mag[:,i]/cs_areas[:,i]
#        cs_surf_hvel_mag[:,i] = cs_surf_hvel_mag[:,i]/cs_areas[:,i]
#        hvel_flux[:,i] = davg_hvel_dir[:,i]   
#        davg_hvel_dir[:,i] = davg_hvel_dir[:,i]/cs_areas[:,i]        
#        cs_surf_hvel_dir[:,i] = cs_surf_hvel_dir[:,i]/cs_areas[:,i]
#
#      ret["davg_salt"] = davg_salt
#      ret["cs_surf_salt"] = cs_surf_salt
#      ret["cs_bott_salt"] = cs_bott_salt
#      ret["davg_den"] = davg_den
#      ret["cs_surf_den"] = cs_surf_den
#      ret["cs_bott_den"] = cs_bott_den
#      ret["davg_temp"] = davg_temp
#      ret["cs_surf_temp"] = cs_surf_temp
#      ret["cs_bott_temp"] = cs_bott_temp
#      ret["davg_hvel_mag"] = davg_hvel_mag
#      ret["cs_surf_hvel_mag"] = cs_surf_hvel_mag
#      ret["hvel_flux"] = hvel_flux
#      ret["davg_hvel_dir"] = davg_hvel_dir
#      ret["cs_surf_hvel_dir"] = cs_surf_hvel_dir
#      ret["names"] = names      
#      ret["dt"] = mr.dt
#      ret["cs_areas"] = cs_areas
#      
#    return ret

  def gen_volume_data(self, base_dir, start_time, end_time, db, data_type="base", hindcast=1):
    ret={}
    
    # use the grid and the shape to get every element in the volume
    if (hindcast):
      [year, week, day, tstep] = model_util.corie2ywdn(start_time)      
      dname = "%s/%d-%02d-%02d/run/" % (base_dir, year, week, db)
      nddir = 7
    else:     #forecast
      [year, day, tstep] = model_util.corie2ydn(start_time)
      dname = "%s/%d-%03d/run/" % (base_dir, year, day)
      nddir = 1
    grid = Gr()
    grid.readHGrid(dname+"hgrid.gr3")

    #collect the elements in the horizontal areas
    indx = where(self.shapes_types == "volume")[0]

    elems = []
    names = []

    for i in indx:
      names.append(self.names[i])
      curr_shape = self.shapes[i]
      elems.append(self.elems_in_area(grid, curr_shape[:,0:2]))

    # get the point data.  

    mr = self.get_model_reader(base_dir, array([]), array([]), start_time, end_time, 0, db, nddir, hindcast)

#    if (self.ob_file):
#      mr = mreader_test(bd = base_dir, bp = array([xpts, ypts]).T, start=start_time, end=end_time, dvg=0, db=db, hindcast=hindcast, self.ob_file, 1)
#    else:
#      mr = mreader_test(bd = base_dir, bp = array([xpts, ypts]).T, start=start_time, end=end_time, dvg=0, db=db, hindcast=hindcast)
    
#    if (data_type=="rho"):
    (sdat, eta, dry) = mr.get_raw_model_data(dtype="salt.63")
    (tdat, eta, dry) = mr.get_raw_model_data(dtype="temp.63")
    depths = mr.grid.depth

    sareas=[]
    # get the element areas
    for i in range(len(elems)):
      tareas = zeros(elems[i].shape[0])
      for j in range(elems[i].shape[0]):
        tareas[j] = signa(mr.grid.x[elems[i][j, 0]], mr.grid.x[elems[i][j, 1]], mr.grid.x[elems[i][j, 2]], \
                          mr.grid.y[elems[i][j, 0]], mr.grid.y[elems[i][j, 1]], mr.grid.y[elems[i][j, 2]])
      sareas.append(tareas)

    # assemble summarized data
    vol_int_salt = zeros((len(elems), eta.shape[1]), float)
    vol_avg_salt = zeros((len(elems), eta.shape[1]), float)
    vol_int_temp = zeros((len(elems), eta.shape[1]), float)
    vol_avg_temp = zeros((len(elems), eta.shape[1]), float)    
    volume = zeros((len(elems), eta.shape[1]), float)

    for i in range(eta.shape[1]):
      z = mr.computeZlevels(depths, eta[:,i])
      dasalt = mr.depth_avg(sdat[:,:,0,i], z)
      datemp = mr.depth_avg(tdat[:,:,0,i], z)      
      for j in range(len(elems)):
        for k in range(elems[j].shape[0]):
          nodes = [elems[j][k, 0], elems[j][k, 1], elems[j][k, 2]]
          if (dry[nodes[0],i]==0 and dry[nodes[1],i]==0 and dry[nodes[2],i]==0):
            dpth = array([(depths[nodes[0]]+eta[nodes[0],i]), (depths[nodes[1]]+eta[nodes[1],i]), (depths[nodes[2]]+eta[nodes[2],i])], float)
            volume[j, i] += sum(dpth)/3*sareas[j][k]

            vol_int_salt[j,i] += sum([dasalt[nodes]*dpth])/3*sareas[j][k]
            vol_int_temp[j,i] += sum([datemp[nodes]*dpth])/3*sareas[j][k]            
        vol_avg_salt[j,i] = vol_int_salt[j,i]/volume[j,i]
        vol_avg_temp[j,i] = vol_int_temp[j,i]/volume[j,i]        

    ret["vol_avg_salt"] = vol_avg_salt
    ret["vol_int_salt"] = vol_int_salt
    ret["vol_avg_temp"] = vol_avg_temp
    ret["vol_int_temp"] = vol_int_temp    
    ret["volume"] = volume
    ret["names"] = names      
    ret["dt"] = mr.dt

    return ret

  def gen_channel_norms(self, base_dir, start_time, end_time, shape_type, db, ndays_dir=7, hindcast=1):
#    if shape_type=="cross_sec" or shape_type=="transect":
    indx = where(self.shapes_types == shape_type)[0]
    for i in indx:
      curr_shape = self.shapes[i]
      mr = self.get_model_reader(base_dir, curr_shape[:,0], curr_shape[:,1], start_time, end_time, 0, db, ndays_dir, hindcast, 1)      
      
      (vdat, eta, dry) = mr.get_model_data(dtype="hvel.64")
      depths = dot(mr.ob.H, mr.grid.depth)      

      davg_u=zeros(eta.shape, float)
      davg_v=zeros(eta.shape, float)
      for j in range(eta.shape[1]):
        z = mr.computeZlevels(depths, eta[:,j])
        davg_u[:,j] = mr.depth_avg(vdat[:,:,0,j], z)
        davg_v[:,j] = mr.depth_avg(vdat[:,:,1,j], z)

      self.shapes[i]= concatenate((self.shapes[i][:,0:2], model_util.norms_from_uv(davg_u, davg_v)), 1)

  def elems_in_area(self, grid, area):
    elems = array([], int)
    for i in range(grid.ne):
      if (pt_inside_poly(array([grid.x[grid.elem[i,2]-1], grid.y[grid.elem[i,2]-1]], float), area) and \
          pt_inside_poly(array([grid.x[grid.elem[i,3]-1], grid.y[grid.elem[i,3]-1]], float), area) and \
          pt_inside_poly(array([grid.x[grid.elem[i,4]-1], grid.y[grid.elem[i,4]-1]], float), area)):
        if (elems.shape[0]==0):
          elems = array(grid.elem[i,2:5], int)
        else:
          elems = vstack((elems, array(grid.elem[i,2:5], int)-1))
    
    return elems

  def save_cs_data(self, data, fname):
    # data is a dictionary object (hash) with data_name : array data pairs
    # stores the data transposed, so that time is rows (making it run up/down the page) and
    #  each cross section is a row
    fid = open(fname, 'w')

    fid.write(str(len(data.keys()))+"\n")
    fid.write(str(data["start_time"]) + " : " + str(data["end_time"]) + "\n") 
    fid.write(str(data["dt"]) + "\n")

    for dtype in data.keys():
      fid.write(dtype+" : ")
      dat = data[dtype]

      if (type(dat)==ndarray):
        if dat.ndim>1:
          fid.write("array : " + str(dat.shape[0])+"\t"+str(dat.shape[1]) +"\n")
          dat = dat.T                
          for i in range(dat.shape[0]):
            for j in range(dat.shape[1]):
              fid.write(str(dat[i,j]))
              if (j < dat.shape[1]-1):
                fid.write("\t")
              else:
                fid.write("\n")
        else:
          fid.write("array :  1\t" + str(dat.shape[0])+"\n")
          for j in range(dat.shape[0]):
            fid.write(str(dat[j])+"\n")
      elif (type(dat)==float):
        fid.write("float\n")
        fid.write(str(dat)+"\n")

      elif (type(dat)==int):
        fid.write("int\n")
        fid.write(str(dat)+"\n")

      elif (type(dat)==list):
        fid.write("list : "+str(len(dat))+"\n")
        for i in range(len(dat)):
          fid.write(str(dat[i]))
          if (i < len(dat)-1):
            fid.write("\t")
          else:
            fid.write("\n")

      elif (type(dat)==str):
        fid.write("str\n"+dat+"\n")

  def read_cs_data(self, fname):
    # data is a dictionary object (hash) with data_name : array data pairs
    # stores the data transposed, so that time is rows (making it run up/down the page) and
    #  each cross section is a row
    ret = {}
    fid = open(fname, 'r')

    ndat = int(fid.readline())
    line = fid.readline().rsplit(':')
    ret["start_time"] = float(line[0])
    ret["end_time"] = float(line[1])
    line = fid.readline()
    ret["dt"] = float(line)

    for i in range(ndat):
      line = fid.readline().rsplit(':')
      name = line[0].strip()

      dtype = line[1].strip()
      if dtype == "array":
        vars = line[2].rsplit('\t')
        ncs = int(vars[0])
        nts = int(vars[1])
        ret[name] = self.read_array_txt(fid, nts, ncs, "float").T
      elif dtype == "list":                         #list of strings 
        vars = fid.readline().rsplit('\t')
        data = []
        for i in range(len(vars)):
          data.append(vars[i].split())        
        ret[name] = data
      elif dtype == "str":
        ret[name] = fid.readline().strip()
      else:
        ret[name] = float(fid.readline().strip())   #currently no "int" or "str" data, although there is a list of str

    return ret

  def read_array_txt(self, fid, nrows, ncols, dtype):
    ret = nan*ones((nrows, ncols), dtype)
    for i in range(nrows):
      line = fid.readline()

      data = line.rsplit("\t")
      if (ncols > 1):
        for j in range(len(data)):
          if dtype == "float":
            ret[i,j] = float(data[j])
          elif dtype == "int":
            ret[i,j] = int(data[j])
          else:
            ret[i,j] = data[j].strip()
      else:
        if dtype == "float":
          ret[i] = float(data[0])
        elif dtype == "int":
          ret[i] = int(data[0])
        else:
          ret[i] = data[0].strip()        

    return ret



########################### scratch ################################

#  gen_cross section norms to return the norms rather than set shapes norms
#    ret = []
#    
#    xpts = array([])
#    ypts = array([])
#
#    if shape_type=="cross_sec" or shape_type=="transect":
#      indx = where(self.shapes_types == "cross_sec")[0]
#
#      xnorms = array([])
#      ynorms = array([])
#      pts_indx = array([], int)
#      for i in indx:
#        names.append(self.names[i])
#        if pts_indx.shape[0]==0:
#          pts_indx = append(pts_indx, 0)
#        else:
#          pts_indx = append(pts_indx, pts_indx[-1]+curr_shape.shape[0])
#        curr_shape = self.shapes[i]
#  
#        xpts = append(xpts, curr_shape[:,0])
#        ypts = append(ypts, curr_shape[:,1])        
#      pts_indx = append(pts_indx, pts_indx[-1]+curr_shape.shape[0])
#
#      #move data collection and norms_from_uv outside a level once pts and (maybe) volumes are added
#      mr = mreader_test(bd = base_dir, bp = array([xpts, ypts]).T, start=start_time, end=end_time, dvg=0, db=db, hindcast=hindcast)
#      (vdat, eta, dry) = mr.get_model_data(dtype="hvel.64")
#      depths = dot(mr.ob.H, mr.grid.depth)      
#
#      davg_u=zeros(eta.shape, float)
#      davg_v=zeros(eta.shape, float)
#      for i in range(eta.shape[1]):
#        z = mr.computeZlevels(depths, eta[:,i])
#        davg_u[:,i] = mr.depth_avg(vdat[:,:,0,i], z)
#        davg_v[:,i] = mr.depth_avg(vdat[:,:,1,i], z)
#
#      for j in range(pts_indx.shape[0]-1):
#        ret.append(array(model_util.norms_from_uv(davg_u[pts_indx[j]:pts_indx[j+1],:], davg_v[pts_indx[j]:pts_indx[j+1],:]), float))
#    return ret
