#!/usr/bin/python
""" This script generates a hotstart.in file based on NCOM 

The option for using different "versions" of NCOM refers to altering the 
salinity values in the ocean by N%.  Estuary values are not affected.

For example: genHotstart(10)
- Generates hotstart.in using NCOM values scaled by 10%

Acceptable values are 0, 10, 15, 20, and 25. 0 is the default value.
"""
import os
import sys
import subprocess
import datetime

SCRIPT_PATH = '/home/workspace/users/lopezj/scripts/selfe_setup/'

# The different versions of ncom refer to the value of the scaled salinity 
# values in the ocean 10%, 15%, 20%, 25%.  They are not scaled in the estuary.
def genHotstart(startDate, ncom_version, grid, n_tracers):
# Generate hotstart.in assuming ocean is only boundary that is of flux type 4.
# First create date.in file
  try:
    datein = open("date.in", "w")
    try:
      if(grid==15):
        datein.write('1 1 1\n')   # 1 open boundary of type 4
      else:
        datein.write('1 1 2\n')   # 2 open boundary of type 4?
      datein.write('%i %i %i 0\n' % 
        (startDate[0], startDate[1], startDate[2]))    
      datein.write('0\n')       # ireduce
      datein.write('Selfe\n')     # selfe or elcirc
      datein.write('1.e6 33.\n')    # ht(depth), salinitymin 
      datein.write('%s\n' % n_tracers)# number of tracers
    finally:
      datein.close()
  except IOError:
    print("error: unable to open date.in for writing\n")
    raise

# Call pragram to generate hotstart.in from ncom values
# The different versions of ncom refer to the value of the scaled salinity 
# values in the ocean.  They are not scaled in the estuary.
  print '\nGenerating hotstart.in\n'
  try:
    _cmd_args = "%s/bin/readncom8b%d" % (SCRIPT_PATH, ncom_version)
    print _cmd_args
    subprocess.check_call(_cmd_args)
  except:
    print("Error: problem generating hotstart.in")
    raise

if __name__ == "__main__":
# grab start and number of days 
  if len(sys.argv) != 7:
    sys.exit("Usage: %s [begin_date mm-dd-yyyy] [number of days in run] "
         "[grid (15, 22, 26, or 27)] [time_step] [number of tracers] "
         "[ncom_version (0,10,20)]"  % sys.argv[0])
  else:
    start = sys.argv[1]
    nDays = int(sys.argv[2])
    grid  = int(sys.argv[3])
    time_step = int(sys.argv[4])
    n_tracers = int(sys.argv[5])
    ncom_version = int(sys.argv[6])

# deal with start and ending days
  (month, day, year) = start.split('-')
  startDay   = int(day)
  startMonth = int(month)
  startYear  = int(year) 

# Python style dates for script.
  startDate = datetime.date(startYear, startMonth, startDay)
  startDate = startDate.timetuple()

# Call to get temp.th and flux.th because temp.th is needed for hotstart.in
# Just to ensure required files are present for generating hotstart.in
  if os.path.isfile('temp.th'): 
    pass
  else:
    print '\nGenertating flux.th and temp.th\n'
    from genFluxTemp  import genFluxTemp
    genFluxTemp(startDate, nDays, grid, time_step)
# Copy *.gr3 files from canonical versions if necessary
# Grid 15 is hi-res.
# Grid 27 is grid 26 cut at Beaver Army.
# Grid else is hi-res 26. 
  if os.path.isfile('hgrid.gr3'):
    pass
  else:
    print '\nlinking to gr3 files for db%i' % grid
    if(grid == 15):
      ret_code = os.system('ln -s %s/hi_res/* ./' % SCRIPT_PATH)
      if ret_code != 0:
        print("Error: problem creating symlinks of grid files to run directory")
        raise IOError
    elif(grid == 22):
      #os.symlink('%s/db22/*' % SCRIPT_PATH, './')
      ret_code = os.system('ln -s %s/db22/* ./' % SCRIPT_PATH)
      if ret_code != 0:
        print("Error: problem creating symlinks of grid files to run directory")
        raise IOError
    elif(grid == 26):
      #os.symlink('%s/db26/*' % SCRIPT_PATH, './')
      ret_code = os.system('ln -s %s/db26/* ./' % SCRIPT_PATH)
      if ret_code != 0:
        print("Error: problem creating symlinks of grid files to run directory")
        raise IOError
    elif(grid == 27):
      ret_code = os.system('ln -s %s/db26c/* ./' % SCRIPT_PATH)
      if ret_code != 0:
        print("Error: problem creating symlinks of grid files to run directory")
        raise IOError
    else: 
# High resolution version of 26, set up grid stuff then the rest is based on 26.
      #os.symlink('%s/db22/*' % SCRIPT_PATH, './')
      ret_code = os.system('ln -s %s/db26_div4/* ./' % SCRIPT_PATH)
      if ret_code != 0:
        print("Error: problem creating symlinks of grid files to run directory")
        raise IOError
      grid = 26

  genHotstart(startDate, ncom_version, grid, n_tracers)

