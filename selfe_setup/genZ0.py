#!/usr/local/bin/python
""" This script is just a wrapper for readssh2b to generate Z0.out

Works for hi-res(15), 22, 26, 26c(27), and hi-res-26(28)
"""

import sys
import os
import datetime
import shutil
import subprocess
import time

# TODO: Add check to make sure hgrid.ll is in directory as well as hgrid.gr3 

# Constants...
BIN_PATH = "/home/workspace/users/lopezj/bin/"

def genZ0(startDate):

# Generate Z0.out for bctides.in.
# This version of readssh2b uses ssh.in for input instead of date.in.
# They are the same thing.
    dateIn = open("./ssh.in", "w")
    str = time.strftime("%Y %m %d", startDate)
    str = str+" 0 \n"
    dateIn.write(str) 
    dateIn.write("0\n")
    dateIn.write("2 1 2\n")         # Number of boundaries that need Z0 
    dateIn.close()

    try:
        _cmd_args = "%s/readssh2b" % BIN_PATH
        subprocess.check_call(_cmd_args, stdout=subprocess.PIPE)
    except subprocess.CalledProcessError, e: 
        print("Error: problem generationg Z0.out from readssh2b")
        raise
    
if __name__ == '__main__':
# Grab start and number of days 
    if len(sys.argv) != 3:
        sys.exit("Usage: %s [begin_date mm-dd-yyyy] [grid (15, 22, 26, or 27)]" 
                 % sys.argv[0])
    else:
        start = sys.argv[1]
        grid  = int(sys.argv[2])
    
    (month, day, year) = start.split('-')
    startDay   = int(day)
    startMonth = int(month)
    startYear  = int(year) 

# Python style dates for script.
    startDate = datetime.date(startYear, startMonth, startDay)
    startDate = startDate.timetuple()

    genZ0(startDate, grid)
