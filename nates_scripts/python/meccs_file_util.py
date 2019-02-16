from numpy import *
import pdb
import sys

def read_pt_file(fname, do_norms):
  try:
    fid = open(fname, 'r')
  except:
    sys.exit("could not open "+fname)

  names = []
  lines = fid.readlines()
  fid.close()
  npts = int(lines[0].rstrip())
  ret = zeros((npts,4),float)  
  for i in range(npts):
    tmp = lines[i+1].split()
    if do_norms:
      ret[i,:] = [float(tmp[0]), float(tmp[1]), float(tmp[2]), float(tmp[3])]
    else:
      ret[i,0:2] = [float(tmp[0]), float(tmp[1])]
    if (len(tmp)==5):
      names.append(tmp[4])
  return (ret, names)


def write_cs_file(fname, cpt, names, rivers = -1):
  try:
    fid = open(fname, 'w')
  except:
    sys.exit("could not open "+fname)

  fid.write(str(len(cpt))+"\n")

  for i in range(len(cpt)):
    npts = len(cpt[i])
#    fid.write(str(npts)+"\t"+names[i]+"\n")
    fid.write(str(npts)+"\t"+names[i])

    if (rivers == -1):
      fid.write("\n")
    else:
      for j in range(len(rivers[i])):
        fid.write("\t"+str(rivers[i][j]))
      fid.write("\n")
    
    for j in range(npts):
      for k in range(3):
        fid.write(str(cpt[i][j][k])+"\t")
      fid.write(str(cpt[i][j][3])+"\n")

  fid.close()
      
def file_exists(f):
  try:
    file=open(f)
  except IOError:
    exists = 0
  else:
    exists = 1
    file.close()
  return exists

def read_cs_file(fname, do_norms=1):
  try:
    fid = open(fname, 'r')
  except:
    sys.exit("could not open "+fname)

  lines = fid.readlines()
  fid.close()
  ncs = int(lines[0].rstrip())
  ret = []
  names = []
  rivers = []
  curr = 1
  for i in range(ncs):
    tmp = lines[curr].split()
    npts = int(tmp[0])
    if (len(tmp)>1):
      names.append(tmp[1].rstrip())
    if (len(tmp)>2):
      curr_rvrs = []
      for j in range(2,len(tmp)):
        curr_rvrs.append(int(tmp[j]))
      rivers.append(curr_rvrs)
    else:
      rivers.append([-1])
    curr = curr+1
    cpt = zeros((npts,4),float)
    for j in range(npts):
      tmp = lines[curr].split()
      curr = curr+1
      if do_norms:
        cpt[j,:] = [float(tmp[0]), float(tmp[1]), float(tmp[2]), float(tmp[3])]
      else:
        cpt[j,0:2] = [float(tmp[0]), float(tmp[1])]
    ret.append(cpt)

  return (ret, names, rivers)

def read_area_file(fname):
  try:
    fid = open(fname, 'r')
  except:
    sys.exit("could not open "+fname)

  lines = fid.readlines()
  fid.close()
  ncs = int(lines[0].rstrip())
  ret = []
  names = []
  curr = 1
  for i in range(ncs):
    tmp = lines[curr].split()
    npts = int(tmp[0])
    if (len(tmp)>1):
      names.append(tmp[1].rstrip())
    curr = curr+1
    cpt = zeros((npts,2),float)
    for j in range(npts):
      tmp = lines[curr].split()
      curr = curr+1
      cpt[j,:] = [float(tmp[0]), float(tmp[1])]
    ret.append(cpt)
  return ret

def write_raw_data(out_dir, fname, data, names, start_time, time_interval, time_units, data_units):
  fid = open(out_dir+"/"+fname, "w")
  fid.write(fname+"\n")
  fid.write("description:\t")
  for i in range(len(names)):
    fid.write(names[i])
    if (i < len(names)-1):
      fid.write("\t")
    else:
      fid.write("\n")      
  fid.write("start time: "+str(start_time)+"\n")
  fid.write("time interval: "+str(time_interval)+"\n")
  fid.write("time units: "+time_units+"\n")
  fid.write("data units: "+data_units+"\n")
  for i in range(data.shape[0]):
    if (data.ndim > 1):
      fid.write("\t")
      for j in range(data.shape[1]):
        fid.write(str(data[i,j]))
        if (j < data.shape[1]-1):
          fid.write("\t")
        else:
          fid.write("\n")
    else:
      fid.write(str(data[i])+"\n")
  fid.close()

def write_csv(fname, data):
  fid = open(fname, "w")
  
  for i in range(data.shape[0]):
    for j in range(data.shape[1]):
      fid.write(str(data[i,j])) 
      if (j < data.shape[1]-1):
        fid.write(", ")
      else:
        fid.write("\n")               
  fid.close()

def cs_to_csvs(fin, fout_base):
  (pts, names, rivers) = read_cs_file(fin, 1)

  for i in range(len(pts)):
    fid = open(fout_base+"_"+str(i)+".csv", "w")

    for j in range(pts[i].shape[0]):
      for k in range(pts[i].shape[1]):
        fid.write(str(pts[i][j,k]))
        if (k==pts[i].shape[1]-1):
          fid.write("\n")
        else:
          fid.write(", ")

    fid.close()

def read_raw_data(fin):
  #if # of header lines or trailing new-lines changes, will need to adjust ret formation
  ret={}
  fid = open(fin, "r")
  lines = fid.readlines()
  fid.close()
  npts = len(lines)-7

  ret['header'] = lines[0].rstrip()
  ret['description'] = lines[1].split('description:')[1].rstrip()
  ret['start_time'] = lines[2].split('start time:')[1].rstrip()
  ret['time_interval'] = lines[3].split('time interval:')[1].rstrip()
  ret['time_units'] = lines[4].split('time units:')[1].rstrip()
  ret['data_units'] = lines[5].split('data units:')[1].rstrip()    

  array_data = []
  for i in range(6, len(lines)):
    line = lines[i].strip()
    str_data = line.split('\t')
    data = []
    for j in range(len(str_data)):
      data.append(float(str_data[j]))
    array_data.append(data)

  ret['data'] = array(array_data)
  return ret

def col_stack_data(data_array):
  #assumes the data is in from the "read_raw_data" "write_raw_data" format
  ret = {}

  for i in range(len(data_array)):
    if i==0:
      data = data_array[i]['data']
      header = data_array[i]['header']
      description = data_array[i]['description']
      start_time = data_array[i]['start_time']
      time_interval = data_array[i]['time_interval']
      time_units = data_array[i]['time_units']
      data_units = data_array[i]['data_units']
    else:
      if (data.shape[0] != data_array[i]['data'].shape[0]):
        print "warning: incompatible data array sizes\n"

      data = column_stack((data, data_array[i]['data']))

      header = header + "\t" + data_array[i]['header']
      description = description + "\t" + data_array[i]['description']

#      if (start_time != data_array[i]['start_time']):
#        print "warning: data start times:"+start_time+", "+data_array[i]['start_time']+"not equal in col_stack_data\n"
      if (time_interval != data_array[i]['time_interval']):
        print "warning: data time_interval:"+time_interval+", "+data_array[i]['time_interval']+"not equal in col_stack_data\n"
      if (time_units != data_array[i]['time_units']):
        print "warning: data time units:"+time_units+", "+data_array[i]['time_units']+"not equal in col_stack_data\n"
      if (data_units != data_array[i]['data_units']):
        print "warning: data units:"+data_units+", "+data_array[i]['data_units']+"not equal in col_stack_data\n"    

  ret['header'] = header
  ret['description'] = description
  ret['start_time'] = start_time
  ret['time_interval'] = time_interval
  ret['time_units'] = time_units
  ret['data_units'] = data_units
  ret['data'] = data

  return ret
