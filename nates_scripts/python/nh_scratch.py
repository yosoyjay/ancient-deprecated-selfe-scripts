  def gen_cross_section_data(self, base_dir, start_time, end_time, db, data_type="base", do_norms=0, ndays_dir=7, hindcast=1):
    ret = {}

    #collect the points in the cross sections:
    indx = where(self.shapes_types == "cross_sec")[0]
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
      
      if (do_norms):
        xnorms = append(xnorms, curr_shape[:,2])
        ynorms = append(ynorms, curr_shape[:,3])        
    pts_indx = append(pts_indx, pts_indx[-1]+curr_shape.shape[0])

    if data_type == "base":     #a basic package of: cs area, depth/surface/bottom avg salt,
                                #hvel (along-channel and magnitude), temp, density

      mr = mreader(bd = base_dir, bp = array([xpts, ypts]).T, start=start_time, end=end_time, dvg=0, db=db, hindcast=hindcast)
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
      
      depths_ts = zeros(eta.shape, float)
      for i in range(eta.shape[1]):
        depths_ts[:,i] = depths + eta[:,i]
        z = mr.computeZlevels(depths, eta[:,i])

        #deal with dry nodes:
        idry = where(dry[:,i] == 1)[0]
        depths_ts[idry, i] = 0

        #get the cross section data
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
              temp1 = sum((tdat[k, bindx[k]:-1, 0, i] + tdat[k, bindx[k]+1:, 0, i])/2)*vSegLens1)
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
      
    return ret
