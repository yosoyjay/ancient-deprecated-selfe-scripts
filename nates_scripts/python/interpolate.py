from numpy import *
import pdb

small = 0.00001

def interp_z_abs(data, depths, z):
  #data is depth values in rows, timesteps in columns
  #depths is absolute depths for data extraction in the same units as z
  #z is depths for data, with the same dimensions as data

  ret = zeros((depths.shape[0], data.shape[1]), float)
  iratios = array([0.,0.])

  for i in range(data.shape[1]):
    for j in range(depths.shape[0]):
      ref_d = z[z.shape[0]-1, i] - depths[j]
      for k in range(data.shape[0]-1):
        if (z[k,i]<=ref_d and z[k+1,i]>ref_d) or (depths[j]==0 and i+2==z.shape[0]):
          iratios[0] = 1-(depths[j]-z[k, i]) / (z[k+1,i]-z[k, i])
          iratios[1] = 1-(z[k+1,i]-depths[j]) / (z[k+1,i]-z[k, i])
          ret[j,i] = data[k,i]*iratios[0]+data[k+1,i]*iratios[1]
          break

  return ret

def interp_z_ratio(data, ratios, depth, zmsl):
  #data is depth values in rows, timesteps in columns
  #ratios is array of ratios 0-1, organized from bottom (1) to top (0)
  #depth is the depth of the water column at msl
  #zmsl is depths at msl, numbered deepest at index 0, surface at last index

#  ratios = sort(ratios.sort()
#  ratios = ratios.fliplr()
  curr_r = 0

  rdepths = depth*ratios

  indexs=zeros((ratios.shape[0], 2), int)
  iratios=zeros((ratios.shape[0], 2), float)

  for i in range(zmsl.shape[0]-1):
    if (zmsl[i]<=rdepths[curr_r] and zmsl[i+1]>rdepths[curr_r]) or (ratios[curr_r]==0 and i+2==zmsl.shape[0]):
      indexs[curr_r, 0] = i
      indexs[curr_r, 1] = i+1
      iratios[curr_r, 0] = 1-(rdepths[curr_r]-zmsl[i]) / (zmsl[i+1]-zmsl[i])
      iratios[curr_r, 1] = 1-(zmsl[i+1]-rdepths[curr_r]) / (zmsl[i+1]-zmsl[i])

      curr_r = curr_r+1
      if (curr_r >= rdepths.shape[0]):
        break

  ret = zeros((ratios.shape[0], data.shape[1]), float)
  for i in range(data.shape[1]):
#    pdb.set_trace()
    ret[:,i] = data[indexs[:,0],i]*iratios[:,0] + data[indexs[:,1],i]*iratios[:,1] 
  
  return ret
