# classes for reading SELFE/ELCIRC data
#
# contains classes Gr and sz_header
# based on elio matlab by Sergey Frolov
#   adapted to python by Nate Hyde

from numpy import *
import struct
import pdb
import geometry
from array import array as _array_
import sys

class Gr:
  '''a class for reading and storing elcirc & selfe grids'''
  def __init__(self):
    self.fname = ""
    self.name = []
    self.np = []
    self.ne = []
    self.nodes = []
    self.elem = []
    self.nn = []
    self.x = []
    self.y = []
    self.depth = []    

  def readHGrid(self, fname):
    file = open(fname, 'r')
    lines = file.readlines()

    self.fname = fname
    self.readHeader(lines)
    self.readNodes(lines)
    if self.type == "gb":
      self.readElem(lines)
      self.readEofLines(lines)
    else:
      self.readEofLines(lines)
    self.flag = 1

  def readHeader(self, lines):
    self.name=lines[0].rstrip()
    tmp = lines[1].split()
    if len(tmp) == 1:
      np = int(tmp[0])
      ne = 0
      type = 'bp'     #build points
    elif len(tmp) == 2:
      ne = int(tmp[0])
      np = int(tmp[1])
      type = 'gb'     #grid with possibly a boundary info
    else:
      print 'error reading ne np'
    self.type = type
    self.ne   = ne
    self.np   = np


  def readNodes(self, lines):
    np = self.np
    nodes = ones([np,4])*nan
    for i in range(np):
      vals = lines[i+2].split()
      nodes[i,0:4] = [int(vals[0]), float(vals[1]), float(vals[2]), float(vals[3])]
    self.nodes = nodes
    self.nn=nodes[:,0]
    self.x=nodes[:,1]
    self.y=nodes[:,2]
    self.depth=nodes[:,3]
    
  def readElem(self, lines):
    ne = self.ne
    np = self.np
    elem = ones([ne,6])*nan
    for i in range(ne):
      line = lines[np+2+i].split()
      if len(line) == 6:
        elem[i,0:6] = [int(line[0]), int(line[1]), int(line[2]), int(line[3]), int(line[4]), int(line[5])]
      elif len(line) == 5:
        elem[i,0:5] = [int(line[0]), int(line[1]), int(line[2]), int(line[3]), int(line[4])]
      else:
        print 'error reading element '+str(i)
    self.elem = elem

  def readEofLines(self, lines):  
    bndStart = 2+self.np+self.ne
    if len(lines) > bndStart:
      self.eofLines = lines[2+self.np+self.ne:len(lines)]
    else:
      self.eofLines   = []
      self.type       = 'gr'     #grid with no boundary info


  def findParentElemsNodesAcos(self, x1, y1):
    ret = zeros((x1.shape[0], 7), float) #return is columns of: [elem #, node#1, node#2, node#3, acos#1, acos#2, acos#3]

    for i in range(self.ne):
      aa = zeros(x1.shape[0], float)
      ar3 = zeros((x1.shape[0], 3), float)

      for j in range(3):
        j1=j+1
        j2=j+2
        if(j1>2):
          j1=j1-3
        if(j2>2):
          j2=j2-3
        in0=int(self.elem[i,j+2])-1
        in1=int(self.elem[i,j1+2])-1
        in2=int(self.elem[i,j2+2])-1
        if j==1:
          ar = float(abs(geometry.signa(self.x[in1],self.x[in2],self.x[in0],self.y[in1],self.y[in2],self.y[in0])))
          if ar<=0:
            print "error: negative area: "+str(ar)+" element "+str(i)

        for k in range(x1.shape[0]):
          ar3[k,j]=geometry.signa(self.x[in1],self.x[in2],x1[k],self.y[in1],self.y[in2],y1[k])
          aa[k]=aa[k]+abs(ar3[k,j])
        
      ae=abs(aa-ar)/ar

      for k in range(x1.shape[0]):
        if ae[k]<=0.00001:
          ret[k,0] = i
          ret[k,1:4] = [self.elem[i,2]-1, self.elem[i,3]-1, self.elem[i,4]-1]          
          ret[k,4:7] = ar3[k,0:3]/ar
          ret[k,4]=max([0,min([1,ret[k,4]])])
          ret[k,5]=max([0,min([1,ret[k,5]])])          
          if(ret[k,4]+ret[k,5]>1):
            ret[k,6]=0
            ret[k,5]=1-ret[k,4]
          else:
            ret[k,6]=1-ret[k,4]-ret[k,5]
    return ret



  def printGridStats(self):
    print "fname: "+self.fname+":\n"
    print "name:  "+self.name+"\n"
    print "np:    "+str(self.np) +"\n"
    print "ne:    "+str(self.ne) +"\n"
    raw_input("Press Enter to continue")
    print "nodes: \n"
    print self.nodes
    print "\n"
    raw_input("Press Enter to continue")
    print self.elem
    print "\n"
    raw_input("Press Enter to continue")
    print "nn:    \n"
    print self.nn
    print "\n"
    print "x:     "+str(len(self.x))+"\n"
    raw_input("Press Enter to continue")    
    print self.x
    print "\n"
    print "y:     "+str(len(self.y))+"\n"
    raw_input("Press Enter to continue")
    print self.y
    print "\n"
    raw_input("Press Enter to continue")    
    print "depth: \n"
    print self.depth
    print "\n"   

class sz_header:
  class _struct:
    pass
  def __init__(self, fname):
    self.fname = fname    
    self.vgrid = self._struct()
    self.hgrid = self._struct()
    fid = open(fname, 'rb')
    self.readHHeader(fid)
    self.readVgrid(fid)
    self.readHgrid(fid)
    self.dataStartPos = fid.tell()
    self.computeStepSize()
    self.computeStepIdx()
    fid.close()

  def readHHeader(self, fid):
    self.dataFormat = str(fid.read(48))
    self.version = str(fid.read(48))
    self.startTime = str(fid.read(48))
    self.varType = str(fid.read(48))
    self.varDimension = str(fid.read(48))
    self.nSteps = struct.unpack('i',fid.read(4))[0]
    self.dt = struct.unpack('f',fid.read(4))[0]
    self.skip = struct.unpack('i',fid.read(4))[0]
    self.flagSv = struct.unpack('i',fid.read(4))[0]
    self.flagDm = struct.unpack('i',fid.read(4))[0]

  def readVgrid(self, fid):
    self.vgrid.startPos = fid.tell()
    self.vgrid.nLevels = struct.unpack('i', fid.read(4))[0]
    self.vgrid.kz = struct.unpack('i', fid.read(4))[0]

    #h0 doesn't seem to be stored, or is not stored in this order.  Hwver, nLevels, startPos, and zLevels and sLevels are right, so
    # the problem may not be major.  Still, something to look into...
    temp = struct.unpack('f', fid.read(4))[0]
    self.vgrid.h0 = 0.1
    self.vgrid.hs = struct.unpack('f', fid.read(4))[0]
    self.vgrid.hc = struct.unpack('f', fid.read(4))[0]
    self.vgrid.theta_b = struct.unpack('f', fid.read(4))[0]
    self.vgrid.theta_f = struct.unpack('f', fid.read(4))[0]
    self.vgrid.zLevels = array([])
    temp =_array_('f')
    temp.fromfile(fid, self.vgrid.kz-1)
    self.vgrid.zLevels = array(temp)
    temp = _array_('f')
    temp.fromfile(fid, self.vgrid.nLevels-self.vgrid.kz+1)
    self.vgrid.sLevels = array(temp)

  def readHgrid(self, fid):
    self.hgrid.startPos = fid.tell();
    self.hgrid.type = 'gr'; 		#it is a grid without a boundary
    self.hgrid.np = struct.unpack('i', fid.read(4))[0]
    self.hgrid.ne = struct.unpack('i', fid.read(4))[0]
    sp = fid.tell()

    hgridTmp1 = _array_('f')
    hgridTmp1.fromfile(fid, 4*self.hgrid.np)
    hgridTmp = array(hgridTmp1)
    hgridTmp = reshape(hgridTmp, (self.hgrid.np,4))
    self.hgrid.x = hgridTmp[:,0];
    self.hgrid.y = hgridTmp[:,1];
    self.hgrid.depth = hgridTmp[:,2];

    fid.seek(sp);
    hgridTmp1 = _array_('i')
    hgridTmp1.fromfile(fid,4*self.hgrid.np)
    hgridTmp = array(hgridTmp1)
    hgridTmp = reshape(hgridTmp, (self.hgrid.np,4))
    self.hgrid.bIdx = hgridTmp[:,3]-1;	#bottom level index
    self.hgrid.nodes = array([range(self.hgrid.np), self.hgrid.x, self.hgrid.y, self.hgrid.depth])
    self.readElem_3(fid)

  def readElem_3(self, fid):
    #assumes all elements are size 3 (i.e. triangles).  Using a mix of shapes requires
    # checking each element size, causing the file read to be very time consuming
    tmp = _array_('i')
    tmp.fromfile(fid, 4*self.hgrid.ne)
    self.hgrid.elem = concatenate((reshape(range(1,self.hgrid.ne+1),(self.hgrid.ne,1)),reshape(tmp,(self.hgrid.ne,4))),1)

  def computeStepSize(self):
    if self.flagDm == 3:
      self.bIdx = maximum(0, self.hgrid.bIdx)
      self.gridSize = sum(self.vgrid.nLevels - (self.bIdx+1) + 1)
    else:
      self.gridSize = self.hgrid.np
    self.stepSize = 2*4 + self.hgrid.np*4 + self.gridSize*4*self.flagSv

  def computeStepIdx(self):
    self.idx = self._struct()
    if self.flagDm == 3:
      #NaN would be better, but Python seems to have limited NaN support
      xv_x = -999999*ones((self.hgrid.np, self.vgrid.nLevels), int) 
      xv_v = -999999*ones((self.hgrid.np, self.vgrid.nLevels), int)
      n_x = arange(self.hgrid.np)
      nvrt = self.vgrid.nLevels

      for i in range(nvrt):
        idx = where(self.hgrid.bIdx <= i)
        xv_x[idx,i] = n_x[array(idx)]
        xv_v[idx,i] = i*ones((array(idx).shape[0],), int)

      idxNodes = reshape(xv_x, (xv_x.shape[0]*xv_x.shape[1],1))
      idxLev = reshape(xv_v, (xv_v.shape[0]*xv_v.shape[1],1))

      idx = where(idxNodes != -999999)
      self.idx.idxNodes = idxNodes[idx]
      idx = idxLev != -999999
      self.idx.idxLev = idxLev[idx]
    if self.flagDm == 2:
      self.idx.idxNodes = arange(self.hgrid.np)
      self.idx.idxLev = array([])

    #precompute mapping from selfe binary to selfe hotstart
    self.idx.idx_all = zeros((self.hgrid.np, self.vgrid.nLevels), int)
    idx_eb = arange(self.idx.idxNodes.shape[0])
    for i in range(self.vgrid.nLevels):
      idxl = where(self.idx.idxLev == i)
      idxn = self.idx.idxNodes[idxl]
      self.idx.idx_all[idxn, i] = idx_eb[idxl]+1 #+1 for logical_and to make self.idx.idx_all_mask:index is one up...

    self.idx.idx_all = reshape(self.idx.idx_all, (multiply.reduce(self.idx.idx_all.shape),), 1)
    self.idx.idx_all_mask = logical_and(self.idx.idx_all, ones(self.idx.idx_all.shape, int))
    self.idx.idx_all[self.idx.idx_all_mask] = self.idx.idx_all[self.idx.idx_all_mask]-1
      
  def readTimeStep(self, tstep):
    if tstep > self.nSteps:
      sys.exit("time step out of range: "+str(tstep)+" > "+str(self.nSteps))

    fid = open(self.fname, 'rb')
      
    data = zeros((self.gridSize, self.flagSv), float)

    fid.seek(self.dataStartPos + self.stepSize*(tstep-1))
    ts = self._struct()
    temp = _array_('f')

    temp.fromfile(fid, 1)
    ts.t = temp[0]
    temp = _array_('i')
    temp.fromfile(fid, 1)
    ts.tit = temp[0]
    temp = _array_('f')
    temp.fromfile(fid, self.hgrid.np)
    ts.eta = array(temp)
    temp = _array_('f')
    temp.fromfile(fid, self.flagSv*self.gridSize)
    data = array(temp)
    if self.flagSv==2:
      data = array([data[::2], data[1::2]]).T
    else:
      data = data.T
    fid.close()

    return (data, ts)

  def map_sz2hts(self, eb_data):
    data_hts = nan*ones(self.idx.idx_all.shape)
    data_hts[self.idx.idx_all_mask]=eb_data[self.idx.idx_all[self.idx.idx_all_mask]]
    data_hts = reshape(data_hts, [self.hgrid.np, self.vgrid.nLevels], 1)
    return data_hts
   
#  def printHeader(self):
#    print "dataFormat: "+self.dataFormat
#    print "\n"
#    print "version: "+self.version
#    print "\n"    
#    print "startTime: "+self.startTime
#    print "\n"    
#    print "varType: "+self.varType
#    print "\n"
#    print "varDimension: "+self.varDimension
#    print "\n"    
#    print "nSteps: "+str(self.nSteps)
#    print "\n"    
#    print "dt: "+str(self.dt)
#    print "\n"    
#    print "skip: "+str(self.skip)
#    print "\n"    
#    print "flagSv: "+str(self.flagSv)
#    print "\n"
#    print "flagDm: "+str(self.flagDm)
#
#  def printVgrid(self):
#    print "startPos: "+str(self.vgrid.startPos)
#    print "nLevels: "+str(self.vgrid.nLevels)
#    print "kz: "+str(self.vgrid.kz)
#    print "h0: "+str(self.vgrid.h0)
#    print "hs: "+str(self.vgrid.hs)
#    print "hc: "+str(self.vgrid.hc)
#    print "theta_b: "+str(self.vgrid.theta_b)
#    print "theta_f: "+str(self.vgrid.theta_f)
#    print "zLevels: "+str(self.vgrid.zLevels)
#    print "sLevels: "+str(self.vgrid.sLevels)         
