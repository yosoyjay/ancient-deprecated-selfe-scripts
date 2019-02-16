#class for creating a matrix used for computing a collection of point values from grid data
#generally ob points are npts X nnodes, and grids are nnodes X nlevels

# ob is based on the elio matlab code by Sergey Frolov
#   adapted to python by Nate Hyde

from numpy import *
import struct
import pdb
import model_util
import gr
from array import array as _array_
import struct

class Ob:
  def __init__(self):
    pass

  def file_init(self, fname, grid):
    self.grid = grid
    self.load_ob(fname)

  def grid_init(self, grid, xy):
    self.grid=grid   # a gr object
    self.x = xy[:,0]
    self.y = xy[:,1]
    self.comp_wxy()
    self.h_indx=arange(self.grid.x.shape[0])    
    self.obxy2h()

  def comp_wxy(self):
    px = self.x
    py = self.y
    pnp = px.shape[0]

    vals = self.grid.findParentElemsNodesAcos(px, py)
    self.tridx = array(vals[:,0]).astype(int)
    self.pidx = array(vals[:,1:4]).astype(int)
    self.w = vals[:,4:7] 

  def obxy2h(self):
    self.H = zeros((self.x.shape[0], self.grid.np), float)
    #could do this more efficiently by indexing arrays with other arrays
    for r in range(self.pidx.shape[0]):
        for c in [0,1,2]:
            self.H[r, self.pidx[r,c]] = self.w[r, c]

    hindx = []
    
    for i in range(self.H.shape[1]):
      if (nonzero(self.H[:,i])[0].shape[0]>0):

        if (len(hindx)==0):
          self.H2 = self.H[:,i]
        else:
          self.H2 = column_stack((self.H2, self.H[:,i]))

        hindx.append(i)
    self.h_indx = hindx

  def ob_h_multiply(self, m):
    if (ndim(m)==1):        
      return model_util.small_dot(self.H2, m[self.h_indx])
    else:        
      return model_util.small_dot(self.H2, m[self.h_indx,:])

  def save_ob(self, fname):
    fid = open(fname, 'wb')
    
    data = struct.pack('i', self.x.shape[0])
    fid.write(data)
    tmp_a = _array_('f', self.x)
    tmp_a.tofile(fid)
    tmp_a = _array_('f', self.y)
    tmp_a.tofile(fid)
    
    tmp_a = _array_('i', self.tridx)
    tmp_a.tofile(fid)
    tmp_a = _array_('i', self.pidx[:,0])
    tmp_a.tofile(fid)
    tmp_a = _array_('i', self.pidx[:,1])
    tmp_a.tofile(fid)
    tmp_a = _array_('i', self.pidx[:,2])
    tmp_a.tofile(fid)
    
    tmp_a = _array_('f', self.w[:,0])
    tmp_a.tofile(fid)
    tmp_a = _array_('f', self.w[:,1])
    tmp_a.tofile(fid)
    tmp_a = _array_('f', self.w[:,2])
    tmp_a.tofile(fid)            

    data = struct.pack('i', self.H.shape[1])
    fid.write(data)        
    for i in range(self.H.shape[0]):
      tmp_a = _array_('f', self.H[i,:])
      tmp_a.tofile(fid)

    data = struct.pack('i', self.H2.shape[1])
    fid.write(data)        
    for i in range(self.H2.shape[0]):
      tmp_a = _array_('f', self.H2[i,:])
      tmp_a.tofile(fid)

    tmp_a = _array_('i', self.h_indx)
    tmp_a.tofile(fid)

  def load_ob(self, fname):
    fid = open(fname, 'rb')
    npts = struct.unpack('i',fid.read(4))[0]

    tmp_a = _array_('f')
    tmp_a.fromfile(fid, npts)
    self.x = array(tmp_a, float)
    tmp_a = _array_('f')
    tmp_a.fromfile(fid, npts)
    self.y = array(tmp_a, float)
    tmp_a = _array_('i')
    tmp_a.fromfile(fid, npts)
    self.tridx = array(tmp_a, int)
    self.pidx = zeros((npts, 3), int)
    tmp_a = _array_('i')
    tmp_a.fromfile(fid, npts)
    self.pidx[:,0] = array(tmp_a, int)
    tmp_a = _array_('i')
    tmp_a.fromfile(fid, npts)
    self.pidx[:,1] = array(tmp_a, int)        
    tmp_a = _array_('i')
    tmp_a.fromfile(fid, npts)
    self.pidx[:,2] = array(tmp_a, int)    

    self.w = zeros((npts, 3), float)
    tmp_a = _array_('f')
    tmp_a.fromfile(fid, npts)
    self.w[:,0] = array(tmp_a, int)
    tmp_a = _array_('f')
    tmp_a.fromfile(fid, npts)
    self.w[:,1] = array(tmp_a, int)        
    tmp_a = _array_('f')
    tmp_a.fromfile(fid, npts)
    self.w[:,2] = array(tmp_a, int)    

    nnodes = struct.unpack('i',fid.read(4))[0]

    self.H = zeros((npts, nnodes), float)
    for i in range(npts):
      tmp_a = _array_('f')
      tmp_a.fromfile(fid, nnodes)
      self.H[i,:] = array(tmp_a)

    ncols = struct.unpack('i',fid.read(4))[0]
    self.H2 = zeros((npts, ncols), float)    
    for i in range(self.H2.shape[0]):
      tmp_a = _array_('f')
      tmp_a.fromfile(fid, ncols)
      self.H2[i,:] = array(tmp_a)

    tmp_a = _array_('i')
    tmp_a.fromfile(fid, ncols)
    self.h_indx = array(tmp_a, int)    


    
    
