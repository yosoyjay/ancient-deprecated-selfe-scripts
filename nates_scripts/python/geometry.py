# basic geometry library

import pdb
import math
from numpy import *

def magnitude(u, v):
  return (u**2 + v**2)**(0.5)

def distance(x1, y1, x2, y2):
  return magnitude(x1-x2, y1-y2)

def pt_inside_poly(pt, poly):
  #if a line extending infinitely in any one direction (2D) crosses the polygon
  # an odd # of times, the point is inside: else outside
  #pt is a single [x,y] point
  #poly is the points in the polygon (assumes the first and last are connected)
  
  pt2 = array([max(poly[:,0])+100, max(poly[:,1])+100], float)

#  pdb.set_trace()

  nintersects = 0
  for i in range(poly.shape[0]-1):
    [int_x, int_y, d] = seg_seg_intersect(array([pt[0], pt2[0]], float), array([pt[1], pt2[1]], float), \
                                          array([poly[i,0], poly[i+1,0]], float), array([poly[i,1], poly[i+1,1]], float))
    if d > 0:
      nintersects = nintersects + 1

  [int_x, int_y, d] = seg_seg_intersect(array([pt[0], pt2[0]], float), array([pt[1], pt2[1]], float), \
                                        array([poly[0,0], poly[-1,0]], float), array([poly[0,1], poly[-1,1]], float))
  if d > 0:
    nintersects = nintersects + 1

  insd = 1
  if nintersects % 2 == 0:
    insd = 0

  return insd

def seg_seg_intersect(lx, ly, sx, sy):

  [a1, b1] = pts2line(lx, ly)
  [a2, b2] = pts2line(sx, sy)

#  pdb.set_trace()

  [x, y] = lines_intersect(array([a1, a2]), array([b1, b2]))
  if isnan(x):
    d = -1
  else:
    sxs = sort(sx)
    sys = sort(sy)
    lxs = sort(lx)
    lys = sort(ly)
    if (x>=sxs[0] and x<=sxs[1] and y>=sys[0] and y<=sys[1]) and (x>=lxs[0] and x<=lxs[1] and y>=lys[0] and y<=lys[1]):
      d = segment_point_dist_2d(array([sxs[0], sys[0]], float), array([sxs[1], sys[1]], float), array([lxs[0], lys[0]], float))
    else:
      d = -1
      
  return array([x, y, d])

def lines_intersect(a, b):
  #a, b are length == 2 arrays of a and by from line equation y=ax+b
  #a[1 or 2] == nan means the line is vertical
  #x, y are the intersect points
  
  if a[0]==a[1]: #same line or parallel
    x = nan
    y = nan
  elif (isnan(a[0])):
    x = b[0]
    y = x*a[1]+b[1]
  elif (isnan(a[1])):
    x = b[1]
    y = x*a[0]+b[0]
  else:
    x = (b[1]-b[0])/(a[0]-a[1])
    y = a[0]*x+b[0]

#  pdb.set_trace()
    
  return [x, y]

def segment_point_dist_2d(p1, p2, p):
  # not yet implemented! 
  return 1

def pts2line(x, y):
  # converts line defined in length 2 arrays (one each for x pts and y pts) into a and b from y=ax+b
  if (x[0]==x[1]):
    a = nan
    b = x[0]
  else:
    a=(y[0]-y[1])/(x[0]-x[1])
    b=y[0]-x[0]*a

  return [a, b]

def signa(x1, x2, x3, y1, y2, y3):    # area formed by 3 points
  return ((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2
