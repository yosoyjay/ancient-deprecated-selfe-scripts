#!/usr/local/bin/python
"""
This script generates the tidal amplitudes and phases for bctides
based loosely upon Joe Cho's guide.  

It's more or less hardcoded, kludgy, and crappy.

Upated to be less crappy and use an early version of grid objects
lopezj - 8/16/2012
"""

#-------------------------------------------------------------------------------
# Imports 
#-------------------------------------------------------------------------------
import sys
import os
import datetime
import shutil
import subprocess
import errno
from optparse import OptionParser

sys.path.append('/home/workspace/users/lopezj/scripts/python/model/grid/')
from grid import *

#-------------------------------------------------------------------------------
# Constants 
#-------------------------------------------------------------------------------
HOME_PATH    = "/home/workspace/users/lopezj/"
SCRIPT_PATH  = "/home/workspace/users/lopezj/scripts/selfe_setup/"
BIN_PATH     = "/home/workspace/users/lopezj/bin/selfe_setup/"
TIDES_PATH   = "/home/workspace/users/lopezj/scripts/selfe_setup/tides/"

#-------------------------------------------------------------------------------
# Functions 
#-------------------------------------------------------------------------------
def tdloadGrid(grid):
    """ Loads a grid into a standard object for parsing """
    # Check if this is a standard grid
    if isinstance(grid,int):
        if grid == 22:
            grid = readHGrid("%s/db22/hgrid.gr3" % SCRIPT_PATH)
        elif grid == 26:
            grid = readHGrid("%s/db26/hgrid.gr3" % SCRIPT_PATH)
        elif grid == 27:
            grid = readHGrid("%s/db26c/hgrid.gr3" % SCRIPT_PATH)
        else:
            print "I'm unfamiliar with grid %d, please provide a path instead" % grid
    elif isinstance(grid,str):
    # Try loading a specified grid
        try:
            grid = readHGrid(grid)
        except OSError, e:
            if e.errno == errno.EEXIST:
                print "The file %s does not exist. Please check path." % grid
            else:
                print "Problem loading %s, but it appears to exist.  I dunno." % grid

    return grid


def genTidalData(startDate, grid, llPath):
    """ 
        Generates ap_?.dat files that have amplitude and phase for each ? boundary 

        startDate - datetime start date of the run
        grid - grid object created by readHGrid()
        llPath - Path to hgrid.ll file
    """
    # Load grid
    # TODO: Write this in a general way so it can be used everywhere
    grid = tdloadGrid(grid)

    # Create symbolic link for input to ecp program
    try:
        os.symlink(llPath, "fort.14.ll")
    except OSError, e:
        if e.errno == errno.EEXIST:
            print "File %s already exists. " % llPath
        else:
            print "Something wrong happend linking %s to fort.14.ll " % llPath
            raise
        
    # Run ecp to get *.nos8 file 
    try:
        cmd_args = "%secp" % BIN_PATH
        subprocess.check_call(cmd_args)
    except IOError:
        print("Problem creating fort.14.nos8")
        raise
    
    # Needed for intel_get to get ap_?.dat files
    os.system("ln -fs %sintel_deg ./" % TIDES_PATH)
    os.system("ln -fs %sedpac2xy.gr3  ./" % TIDES_PATH)
    os.system("ln -fs %steanl.tct     ./" % TIDES_PATH)
    os.system("ln -fs %sintel_deg.in ./" % TIDES_PATH)
    os.system("cp fort.14.nos8 nos8")               

    # Loop over open boundaries and only write files for those that are tidally forced
    for ob in range(1,grid.obNodes.nB+1):
        var = "grid.obNodes.ob_%d_type" % ob
        if eval(var) is not 'ocean':
            continue
        try:
            var = "grid.obNodes.ob_%d_nodes" % ob
            of = open("nos8_%i.new" % ob, "w")
            shutil.copyfileobj(open("nos8", "r"), of)
            of.write("%i\n" % grid.obNodes.nB)  # Number of open boundaries
            of.write("%i\n" % grid.obNodes.nB)  # Apparently this doesn't matter
            of.write("%i\n" % len(eval(var)))   # Number of nodes on this boundary
            for node in eval(var):
                of.write("%i\n" % int(node))    # Write each node number on this boundary
            of.close()
        except IOError, e:
            print "Problem writing to nos8_%i.new" % ob

        try:
            args = ["cp", "-f", "nos8_%i.new" % ob, "fort.14.nos8"]
            subprocess.check_call(args)
        except :
            print("Error: Problem copying nos8_i% to fort.14.nos8" % ob)
            raise

        # Run genbc to get *.sta 
        try:
            _cmd_args = ["%sgenbcs" % TIDES_PATH]
            subprocess.check_call(_cmd_args)
        except :
            print("Error: Problem with call to genbcs.f")
            raise

        # Eliminated prompt/response version of intel_deg 
        # TODO: Move from os.system to subprocess.Popen
        os.system("./intel_deg &> /dev/null")
        os.system("mv ap.out ap_%d.dat" % ob)

#-------------------------------------------------------------------------------
# Main 
#-------------------------------------------------------------------------------
if __name__ == '__main__':
    usage = ("%prog [startDate mm-dd-yyyy]) [Path to hgrid.ll] and specify path"
             " to hgrid.gr3 via -p or specify a grid (22,26) via -g")
             
    parser = OptionParser(usage=usage)
    parser.add_option("-p", "--path", action="store", type="string",
                      dest="path", help="Path to hgrid.gr3")
    parser.add_option("-g", "--grid", action="store", type="int",
                      dest="grid", help="Grid number [22, 26, 27(26 cut at BA)")

    parser.set_defaults(path='')
    parser.set_defaults(grid=999)

    (options, args) = parser.parse_args()
    grid = options.grid
    path = options.path
    if (grid == 999 and path == '') or (grid != 999 and path != ''):
        parser.error("Apply either a path to a grid or a grid number")

# grab start and number of days 
    if len(args) != 2:
        parser.error("Incorrect number of arguments")
    else:
        start = args[0]
        llPath = args[1]

# deal with start days
    (month, day, year) = start.split('-')
    startDay   = int(day)
    startMonth = int(month)
    startYear  = int(year) 
    startDate = datetime.datetime(startYear, startMonth, startDay)

    genTidalData(startDate, grid, llPath)
