# class for reading ELCIRC/SELFEmodel data, given a series of points, database (assumes perticular file structure), and start/end times
#
# by Nate Hyde, with computeZlevels borrowed heavily from Sergey Frolov's elio matlab code
#
# Note that Ob initialization takes a good deal of time.  May want to make use of Ob_object.save_ob(fname) and Ob_object.file_init()

from numpy import *
from gr import *
import model_util
from ob import *
#import pdb

class mreader_test:
  def __init__(self, bd="", bp=array([]), start=2300, end=2301, dvg=0, db=16, ndays_dir=7, hindcast=1, ob_file="", do_ob_file=0):
    self.base_dir = bd
    self.base_pts = bp
    self.start_time = start
    self.end_time = end
    self.d_avg = dvg
    self.db = db
    self.bt_indx = []   #index of the surface and bottom for each of the points (self.base_pts)
    self.ndays_per = ndays_dir  #number of model days in each directory
    self.hindcast = hindcast

    if (self.hindcast):
      [year, week, day, tstep] = model_util.corie2ywdn(self.start_time)      
      dname = "%s/%d-%02d-%02d/run/" % (self.base_dir, year, week, self.db)
    else:     #forecast
      [year, day, tstep] = model_util.corie2ydn(self.start_time)
      dname = "%s/%d-%03d/run/" % (self.base_dir, year, day)
    self.grid = Gr()
    self.grid.readHGrid(dname+"hgrid.gr3")

    self.ob = Ob()

    if (bp.shape[0]>0):
      if (do_ob_file==0):
        self.ob.grid_init(self.grid, self.base_pts)
#        self.ob.save_ob(ob_file)
      else:
        self.ob.file_init(ob_file, self.grid)
    
    hdr = sz_header(dname+"1_elev.61")
    self.grid.vgrid = hdr.vgrid
    self.dt = hdr.dt

  def depth_avg(self, data, z):
    dat = (data[:,0:data.shape[1]-1] + data[:,1:data.shape[1]])/2
    depths = z[:,1:z.shape[1]]-z[:,0:z.shape[1]-1]
    dat_x_d = dat*depths
    dat_x_d[isnan(dat_x_d)]=0

    #prevent divide-by-zero errors
    dsums = sum(depths*logical_not(isnan(dat)),1)
    nsums = sum(dat_x_d,1)
    indx = where(dsums==0)
    dsums[indx]=1
    nsums[indx]=0
    
    ret = nsums/dsums

    return ret

  def computeZlevels(self, dp, eta):
    nZ = self.grid.vgrid.zLevels.shape[0]    
    ret = zeros((eta.shape[0], self.grid.vgrid.sLevels.shape[0]+nZ), float)

    #Z levels for Z portion
    for i in range(nZ):
      ret[:,i] = self.grid.vgrid.zLevels[i]

    #Z levels for sigma grid portion
    for i in range(self.grid.vgrid.sLevels.shape[0]):
      sigma = self.grid.vgrid.sLevels[i]      
      C = (1-self.grid.vgrid.theta_b)*sinh(self.grid.vgrid.theta_f*sigma)/sinh(self.grid.vgrid.theta_f) + self.grid.vgrid.theta_b*(tanh(self.grid.vgrid.theta_f*(sigma+0.5))-tanh(self.grid.vgrid.theta_f/2) )/ (2*tanh(self.grid.vgrid.theta_f/2))

      #special case
      idx = where(dp <= self.grid.vgrid.hc)
      ret[idx[0],i+nZ] = sigma*(dp[idx]+eta[idx])+eta[idx]

      #normal case
      idx = where(dp > self.grid.vgrid.hc)      
      ret[idx[0],i+nZ] = eta[idx[0]]*(1+sigma) + self.grid.vgrid.hc*sigma + (dp[idx[0]]-self.grid.vgrid.hc)*C

      #special case 2
      criterion2 = -self.grid.vgrid.hc-(dp-self.grid.vgrid.hc)*self.grid.vgrid.theta_f/sinh(self.grid.vgrid.theta_f)
      idx = where(eta <= criterion2)
      etam  = 0.98*( -self.grid.vgrid.hc-(dp[idx[0]]-self.grid.vgrid.hc)*self.grid.vgrid.theta_f/sinh(self.grid.vgrid.theta_f) );          
      #    warning('S->Z coordinate transformation may not be valid')
      sigma_hat = (ret[idx[0],i+nZ]-etam)/(dp[idx[0]]+etam)
      ret[idx[0],i+nZ] = sigma_hat*(dp[idx[0]]+eta[idx[0]])+eta[idx[0]]

      #special case of dry nodes
      #zoseph assigns z=0 in this case,
      #i (sergey) think it's more consistent to assign depth value in this case

      idx = where((dp + eta)<=self.grid.vgrid.h0)
      #  z[idx,i]      = 0           #joseph
      ret[idx[0],i+nZ] = dp[idx[0]]        #sergey
      
    return ret

  def get_model_data(self, dtype, start=-99999999, end=-99999999):
    if start != -99999999:
      self.start_time = start
    if end != -99999999:
      self.end_time = end

    if (self.hindcast):
      [year, week, day, tstep] = model_util.corie2ywdn(self.start_time)      
      dname = "%s/%d-%02d-%02d/run/" % (self.base_dir, year, week, self.db)
    else:     #forecast
      [year, day, tstep] = model_util.corie2ydn(self.start_time)
      dname = "%s/%d-%03d/run/" % (self.base_dir, year, day)
      
#    [year, week, day, tstep] = model_util.corie2ywdn(self.start_time)
#    dname = "%s/%d-%02d-%02d/run/" % (self.base_dir, year, week, self.db)

    fname = "%d_%s" %  (1, dtype)

    hdr = sz_header(dname+fname)
    step_len = hdr.dt/86400

#    depths = model_util.small_dot(self.ob.H, self.grid.depth)
    depths = self.ob.ob_h_multiply(self.grid.depth)
    z = self.computeZlevels(depths, zeros((self.ob.H.shape[0])))
    nzlvls = z.shape[1]

    tsteps = arange(self.start_time, self.end_time+step_len, step_len)
    if self.d_avg:
      results = zeros((self.ob.x.shape[0], 1, hdr.flagSv, tsteps.shape[0]), float)
    else:
      results = zeros((self.ob.x.shape[0], hdr.vgrid.nLevels, hdr.flagSv, tsteps.shape[0]), float)      
    eta = zeros((self.ob.x.shape[0], tsteps.shape[0]), float)
    dry = zeros((self.ob.x.shape[0], tsteps.shape[0]), int)

    #control variables:
    index = 0
    max_t_tindex = zeros(self.ob.x.shape[0])
    last_day = self.ndays_per + 1

    for d_time in tsteps:
      if (self.hindcast):
        [year, week, day, tstep] = model_util.corie2ywdn(d_time)      
        dname = "%s/%d-%02d-%02d/run/" % (self.base_dir, year, week, self.db)
        fname = "%d_%s" %  (day, dtype)
          
      else:     #forecast
        [year, day, tstep] = model_util.corie2ydn(d_time)
        dname = "%s/%d-%03d/run/" % (self.base_dir, year, day)
        fname = "%d_%s" %  (1, dtype)      

      if day != last_day:        
        try:
          hdr.fid.close()
        except AttributeError:
          pass

        hdr=sz_header(dname+fname)
        last_day=day

      (dd, ts) = hdr.readTimeStep(tstep)
      if hdr.flagSv==2:
        dat = array([hdr.map_sz2hts(dd[:,0]), hdr.map_sz2hts(dd[:,1])])
        if self.d_avg:
#          results[:,0,:,index] = array([self.depth_avg(model_util.small_dot(self.ob.H, dat[0,:,:]), z), self.depth_avg(model_util.small_dot(self.ob.H,dat[1,:,:]), z)]).T  #averages according to MSL: better way is to
          results[:,0,:,index] = array([self.depth_avg(self.ob.ob_h_multiply(dat[0,:,:]), z), self.depth_avg(self.ob.ob_h_multiply(dat[1,:,:]), z)]).T  #averages according to MSL: better way is to
                                                                                #get current depths.
        else:
#          results[:,:,0,index] = model_util.small_dot(self.ob.H, dat[0,:,:])
#          results[:,:,1,index] = model_util.small_dot(self.ob.H, dat[1,:,:])

          results[:,:,0,index] = self.ob.ob_h_multiply(dat[0,:,:])
          results[:,:,1,index] = self.ob.ob_h_multiply(dat[1,:,:])    
      else:
        dat = hdr.map_sz2hts(dd[:])

        #"Illegal Instruction" crash is a problem w/ large arrays on numpy.dot() on 64 bit machines
        # short-term fix is to break the array down to a smaller size.  The method & sizes used here are mostly
        # arbitrary, based on a just a few trials        
        if self.d_avg:
#          results[:,0,0,index] = self.depth_avg(model_util.small_dot(self.ob.H,dat[:,:]), z)
          results[:,0,0,index] = self.depth_avg(self.ob.ob_h_multiply(dat[:,:]), z)          
        else:
          results[:,:,0,index] = self.ob.ob_h_multiply(dat[:,:])

      eta[:,index] = self.ob.ob_h_multiply(ts.eta)
#      if (index==14):
#        pdb.set_trace()      
      dindx = unique(where((self.grid.depth[self.ob.pidx] + ts.eta[self.ob.pidx]) <= self.grid.vgrid.h0)[0])
      dry[dindx, index] = 1

      index=index+1
    return (results, eta, dry)    
   
  def get_raw_model_data(self, dtype, start=-99999999, end=-99999999):
    if start != -99999999:
      self.start_time = start
    if end != -99999999:
      self.end_time = end

    if (self.hindcast):
      [year, week, day, tstep] = model_util.corie2ywdn(self.start_time)      
      dname = "%s/%d-%02d-%02d/run/" % (self.base_dir, year, week, self.db)
    else:     #forecast
      [year, day, tstep] = model_util.corie2ydn(self.start_time)
      dname = "%s/%d-%03d/run/" % (self.base_dir, year, day)      
    fname = "%d_%s" %  (1, dtype)

    hdr = sz_header(dname+fname)
    step_len = hdr.dt/86400

    depths = self.grid.depth
    z = self.computeZlevels(depths, zeros((depths.shape[0])))
    nzlvls = z.shape[1]

    tsteps = arange(self.start_time, self.end_time+step_len, step_len)
    if self.d_avg:
      results = zeros((self.grid.x.shape[0], 1, hdr.flagSv, tsteps.shape[0]), float)
    else:
      results = zeros((self.grid.x.shape[0], hdr.vgrid.nLevels, hdr.flagSv, tsteps.shape[0]), float)      
    eta = zeros((self.grid.x.shape[0], tsteps.shape[0]), float)
    dry = zeros((self.grid.x.shape[0], tsteps.shape[0]), int)

    #control variables:
    index = 0
    max_t_tindex = zeros(self.grid.x.shape[0])
    last_day = self.ndays_per + 1

    for d_time in tsteps:
      if (self.hindcast):
        [year, week, day, tstep] = model_util.corie2ywdn(d_time)      
        dname = "%s/%d-%02d-%02d/run/" % (self.base_dir, year, week, self.db)
        fname = "%d_%s" %  (day, dtype)
          
      else:     #forecast
        [year, day, tstep] = model_util.corie2ydn(d_time)
        dname = "%s/%d-%03d/run/" % (self.base_dir, year, day)
        fname = "%d_%s" %  (1, dtype)      

      if day != last_day:        
        try:
          hdr.fid.close()
        except AttributeError:
          pass

        hdr=sz_header(dname+fname)
        last_day=day

      (dd, ts) = hdr.readTimeStep(tstep)
      if hdr.flagSv==2:
        dat = array([hdr.map_sz2hts(dd[:,0]), hdr.map_sz2hts(dd[:,1])])
        if self.d_avg:
          results[:,0,:,index] = array([self.depth_avg(dat[0,:,:], z), self.depth_avg(dat[1,:,:], z)]).T  #averages according to MSL: better way is to
                                                                                #get current depths.
        else:
          results[:,:,0,index] = dat[0,:,:]
          results[:,:,1,index] = dat[1,:,:]
      else:
        dat = hdr.map_sz2hts(dd[:])

        if self.d_avg:
          results[:,0,0,index] = self.depth_avg(dat[:,:], z)
        else:
          results[:,:,0,index] = dat[:,:]

      eta[:,index] = ts.eta

      dindx = unique(where((self.grid.depth + ts.eta) <= self.grid.vgrid.h0)[0])
      dry[dindx, index] = 1

      index=index+1
    return (results, eta, dry)
    
  def get_river_flux(self, base_dir, start_time, end_time, db, dt, hindcast=1):
    # not thoroughly tested: only confident on whole days (starting at 0:15 on a given day) for forecasts

    
    small = 0.00001
    ntstep_day = int(86400/dt)
    curr_time = start_time

    ret = array([])
    iret = 0

    while curr_time < end_time:  
      if (hindcast):
        [year, week, day, tstep] = model_util.corie2ywdn(curr_time)      
        dname = "%s/%d-%02d-%02d/run/" % (base_dir, year, week, db)
      else:     #forecast
        [year, day, tstep] = model_util.corie2ydn(curr_time)
        dname = "%s/%d-%03d/run/" % (base_dir, year, day)    
  
      [rdt, flx] = self.read_flux_th(dname+"flux.th")
      fstps_mstp = dt/rdt

      if (~ret.any()):
        if flx.ndim>1:
          ret = zeros((ceil((end_time-start_time)*ntstep_day+1), flx.shape[1]), float)
        else:
          ret = zeros((ceil((end_time-start_time)*ntstep_day+1), 1), float)          
      
      if (hindcast):
        flx_start = ntstep_day*(day-1) + fstps_mstp*(tstep-1)
        if end_time <= curr_time+(8-day):
          flx_end = flx_start + (end_time-curr_time) * fstps_mstp*ntstep_day
        else:
          flx_end = flx_start + (ceil(curr_time+small)-curr_time) * fstps_mstp*ntstep_day
      else:
        flx_start = fstps_mstp*(tstep-1)
        if end_time <= curr_time:
          flx_end = flx_start + (end_time-curr_time) * fstps_mstp*ntstep_day
        else:
          flx_end = flx_start + fstps_mstp*ntstep_day

      ncurrsteps = int((flx_end-flx_start)/fstps_mstp)

      for i in range(ncurrsteps):
#        ret[iret+i] = mean(sum(flx[flx_start+i*fstps_mstp:flx_start+(i+1)*fstps_mstp, 0:], 1))
#        ret[iret+i] = mean(sum(flx[flx_start+i*fstps_mstp:flx_start+(i+1)*fstps_mstp, river_nums], 1))
        ret[iret+i,:] = mean(flx[flx_start+i*fstps_mstp:flx_start+(i+1)*fstps_mstp, :], 0)

      iret = iret+ncurrsteps
      curr_time = curr_time + flx_end/(fstps_mstp*ntstep_day)
    return ret
      
  def read_flux_th(self, fname):
    fid = open(fname, 'r')
    
    lines   = fid.readlines()
    nrvr = len(lines[0].split())-1
    ret = zeros((len(lines), nrvr), float)

    dt = int(float(lines[0].split()[0]))

    for i in range(len(lines)):
      vals = lines[i].split()
      for j in range(1, len(vals)):
        ret[i,j-1] = -float(vals[j])

    return [dt, ret]  
