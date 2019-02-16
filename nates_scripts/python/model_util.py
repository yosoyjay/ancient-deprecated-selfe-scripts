# useful functions when working with SELF/ELCIRC data
#
# written by Nate Hyde, borrowing liberally from Sergey Frolov's elio matlab code

import pdb
import math
import geometry
from numpy import *
import sys

def isleap(year):
  return ((year%4==0) and ((year%100!=0) or (year%400==0)))

def days_in_month(month, year):
  days = list([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
  if (isleap(year)):
    days[1] = 29
  return days[month-1]

def days_in_year(year):
  if isleap(year):
    return 366
  else:
    return 365

def ymdhms2corie(yr, mth, d, hr, min, s):
  ret = 0    
  if yr >= 1996:
    for i in range(1996, yr):
      ret=ret+days_in_year(i)
  else:
    pass
#    ret=-1                                #gotta re-think this: later
#    for i in range(yr, 1996):
#      ret=ret-days_in_year(i)   
  for i in range(1,mth):
    ret=ret+days_in_month(i-1, yr)
  ret = ret+d+hr/24+min/(24*60)+s/(24*3600)
  return ret

def ywdhms2corie(yr, wk, d, hr, min, s):   #TODO: negative dates
  ret = 0.
  if yr >=1996:
    for i in range(1996, yr):
      ret=ret+days_in_year(i)
    ret = ret+(wk-1)*7.+d+hr/24.+min/(24.*60.)+s/(24.*3600.)
  return ret

def ydhms2corie(yr, d, hr, min, s):   # TODO: negative dates
  ret = 0.
  if yr >=1996:
    for i in range(1996, yr):
      ret=ret+days_in_year(i)
    ret = ret+d+hr/24.+min/(24.*60.)+s/(24.*3600.)
  return ret

def corie2datevec(dn):          #test this further: may not work on negative dates
  days_count = 0
  year = 1996
  if (dn<1):
    while dn <= days_count-days_in_year(year-1):
      days_count -= days_in_year(year-1)
      year=year-1
    year=year-1
    day = dn - days_count
    month=12
    while abs(day) > days_in_month(month, year):
      day += days_in_month(month, year)
      month -= 1
    day = days_in_month(month, year) + day
  else:
    while int(dn) > days_count+days_in_year(year):
      days_count += days_in_year(year)
      year=year+1
    day = dn - days_count
    month=1

    while int(day) > days_in_month(month, year):
      day -= days_in_month(month, year)
      month += 1
  return [year, month, int(day)]

def corie2ywdn(corie_time):
  dim = list([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
  [year, month, day] = corie2datevec(corie_time)
  part_day = corie_time - math.floor(corie_time)
  doy = sum(dim[0:month-1])+day
  if isleap(year) and month > 2:
    doy=doy+1
  
  week = math.floor((doy-1)/7)+1
  day = (doy-1)%7 + 1
  ts = round(96*part_day)  

  #ts 0 is stored in the previous day
  if ts == 0:
    ts = 96
    day = day-1
    if day==0:
      day = 7
      week = week-1
      if week==0:
        week = 53
        year = year-1
        if isleap(year):
          day = 2
        else:
          day = 1

  return [year, week, day, ts]

def corie2ydn(corie_time):
  dim = list([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])  
  part_day = corie_time - math.floor(corie_time)
  [year, month, day] = corie2datevec(corie_time)  
  doy = sum(dim[0:month-1])+day
  if isleap(year) and month > 2:
    doy=doy+1
  ts = round(96*part_day)  

  #ts 0 is stored in the previous day
  if ts == 0:
    ts = 96
    doy = doy-1
    if doy==0:
      doy = days_in_year(year-1)
      year = year-1

  return [year, doy, ts]

def resolve_line(ptsx, ptsy, vectU, vectV, dist):

  # No error checking!
  # Assumes "dist" is less than or equal the distance between points
  
  xdiff = ptsx[1:]-ptsx[0:-1]
  ydiff = ptsy[1:]-ptsy[0:-1]

  mag = geometry.magnitude(xdiff, ydiff)

  steps = floor(mag/dist)+1

  sz = sum(steps)+1

  ret=zeros((sz, 4), float)

  ret[0,:] = [ptsx[0], ptsy[0], vectU[0], vectV[0]]
  rindx = 1
  for i in range(0,xdiff.shape[0]):
    xdist = xdiff[i]/steps[i]
    ydist = ydiff[i]/steps[i]
    unorm = (vectU[i]+vectU[i+1])/2

    vnorm = (vectV[i]+vectV[i+1])/2

    for j in arange(steps[i]-1):
      ret[rindx, :] = [ret[rindx-1,0]+xdist, ret[rindx-1,1]+ydist, unorm, vnorm];
      rindx = rindx+1

    ret[rindx, :] = [ptsx[i+1], ptsy[i+1], vectU[i+1], vectV[i+1]]
    rindx = rindx+1

  return ret

def density(salt, temp):
  #can give invalid value warnings if an array with (-) values is included
  return 1000-0.157406+0.06793952*temp-0.009095290*temp**2+0.0001001685*temp**3-0.000001120083*temp**4+0.0000000006536332 \
         *temp**5+salt*(0.824493-0.0040899*temp+0.000076438*temp**2-0.00000082467*temp**3+0.0000000053875*temp**4)+ \
         salt**(1.5)*(-0.00572466+0.00010227*temp-0.0000016546*temp**2)+0.00048314*salt**2

def small_dot(A, B, max_ops = nan):
  #small_dot breaks up matrix A to calculate the matrix multiplication dot(A,B)
  # this is necessary because numpy has issues with large matrix multiplication on 64 bit machines
  # simple function only breaks up rows of A

  if isnan(max_ops):
    max_ops = 9000000   # reasonable max size from a few trials on amb25

  nops = A.shape[0]
  if (A.ndim > 1):    
    nops = nops*A.shape[1]
  ncols = 1
  if (B.ndim > 1):
    nops = nops*B.shape[1]
    ncols = B.shape[1]

  if nops <= max_ops:
    return dot(A,B)
  else:
    nsteps = floor(nops / max_ops)
    if (nsteps > A.shape[0]):
      sys.exit("error in model_util.small_dot")
    stepsz = floor(A.shape[0]*(float(max_ops)/float(nops)))
    if (B.ndim > 1):                  #numpy silliness: [npts X 1] size array for B treated as [npts] for assignment
      ret = zeros((A.shape[0], ncols), float)
      for i in range(int(nsteps)):
        ret[i*stepsz:(i+1)*stepsz, :] = dot(A[i*stepsz:(i+1)*stepsz, :], B)
      ret[nsteps*stepsz:nsteps*stepsz+A.shape[0]-nsteps*stepsz, :] = dot(A[nsteps*stepsz:nsteps*stepsz+A.shape[0]-nsteps*stepsz, :], B)
    else:
      ret = zeros((A.shape[0]), float)
      for i in range(int(nsteps)):
        ret[i*stepsz:(i+1)*stepsz] = dot(A[i*stepsz:(i+1)*stepsz, :], B)
      ret[nsteps*stepsz:nsteps*stepsz+A.shape[0]-nsteps*stepsz] = dot(A[nsteps*stepsz:nsteps*stepsz+A.shape[0]-nsteps*stepsz, :], B)        

    return ret
      
def tidal_princ_comp(u,v):
  if (rank(u)==2):
    theta=zeros((u.shape[0]),float)
    var1=zeros((u.shape[0]),float)
    var2=zeros((u.shape[0]),float)
    u1=zeros(u.shape,float)
    v1=zeros(u.shape,float)
    for i in range(u.shape[0]):
      u1[i,:]=u[i,:]-mean(u[i,:])
      v1[i,:]=v[i,:]-mean(v[i,:])      
      mean_u2=mean(u1[i,:]**2)
      mean_v2=mean(v1[i,:]**2)
      mean_uv=mean(u1[i,:]*v1[i,:])
      theta[i]=0.5*arctan2(2*mean_uv, (mean_u2-mean_v2))
      var1[i]=0.5*(mean_u2+mean_v2+((mean_u2-mean_v2)**2 + 4*mean_uv**2)**0.5)
      var2[i]=0.5*(mean_u2+mean_v2-((mean_u2-mean_v2)**2 + 4*mean_uv**2)**0.5)
  else:
    u1=u-mean(u,1)
    v1=v-mean(v,1)
    mean_u2=mean(u1**2)
    mean_v2=mean(v1**2)
    mean_uv=mean(u1*v1)
    theta=0.5*arctan2(2*mean_uv, (mean_u2-mean_v2))
    var1=0.5*(mean_u2+mean_v2+((mean_u2-mean_v2)**2 + 4*mean_uv**2)**0.5)
    var2=0.5*(mean_u2+mean_v2-((mean_u2-mean_v2)**2 + 4*mean_uv**2)**0.5)
  return array([theta, var1, var2], float).T

def norms_from_uv(u1,v1):
  tpc = tidal_princ_comp(u1,v1)

  if (rank(u1)==2):
    ret = array([cos(tpc[:,0]), sin(tpc[:,0])], float).T
  else:
    ret = array([cos(tpc[0]), sin(tpc[0])], float)    
  return ret

def getmonthday(year, dd):
  if ( dd > 366 or dd < 1 ):
    print "Bad yearday = $dd\n";
    return ( -1, -1 );
  if ( year % 4 == 0 and year % 100 ):
    julian_days = [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
  else:
    julian_days = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]

  month = 12;
  while (dd <= julian_days[month]):
    month = month-1
  day = dd - julian_days[month]
  return (month + 1, day);
