#!/usr/local/bin/python
""" Creates bctides.in for a Selfe run

This is flexible to work with different grids, but the open boundaries are 
hard coded and depend on being named in a specific manner in hgrid.gr3. If new
boundaries are added to a grid they should be added to either the RIVER_BOUNDS
or OCEAN_BOUNDS here as well as noting the changes in grid.py which this
code depends on.

"""

#-------------------------------------------------------------------------------
# Imports
#-------------------------------------------------------------------------------
import os
import sys
import glob
import time
import datetime
import linecache
import subprocess
import numpy as np
from optparse import OptionParser

# Hacky --- Need to move things to *.egg  
SCRIPT_PATH = '/home/workspace/users/lopezj/scripts/selfe_setup/'
sys.path.append(SCRIPT_PATH)
from genTidalData import genTidalData, tdloadGrid
from genZ0 import genZ0

sys.path.append('/home/workspace/users/lopezj/scripts/python/model/grid/')
from grid import *

#-------------------------------------------------------------------------------
# Constants 
#-------------------------------------------------------------------------------
# This is NARR, if this fails then we should use NAM = 1
SFLUX_SRC = 3

# Constituent frequencies
N_CONS_FREQS = 9
FREQS ={"Z0":"0.000000e+00",
        "O1":"6.759774e-05",
        "K1":"7.292116e-05",
        "Q1":"6.495457e-05",
        "P1":"7.251056e-05",
        "K2":"1.458423e-04",
        "N2":"1.378797e-04",
        "M2":"1.405189e-04",
        "S2":"1.454441e-04"}

FREQS_NAMES = ['Z0', 'O1', 'K1', 'Q1', 'P1', 'K2', 'N2', 'M2', 'S2']

# River boundary conditions based on precedence.
# This is used to write out the boundaries as a string, but not include
# the number of nodes along the boundary.
# These values are based on what are used for hindcast runs or forecast if no
# hindcast runs are available.
# The names here match those found in grid.py
RIVER_BOUNDS = {'fraser': '0 2 3 3 ! Fraser\n\
0. ! BC discharge\n\
1. ! Relax T\n\
1. ! Relax S\n',
                'beaver': '0 1 1 2 ! Beaver Army\n\
0.\n',
                'bonneville': '0 1 1 2 ! Bonneville Dam\n\
1. ! Relax for T\n\
0. ! BC for salt\n\
1. ! Relax for S\n',
                'morrison':  '0 0 0 0  ! Morrison Bridge\n',
                'willamette': '1 1 1 2 ! Willamette River Falls\n\
1. ! Relax for T\n\
0. ! BC for salt\n\
1. ! Relax for S\n',
                'b_power1': '0 1 1 2 ! Bonneville power1\n\
1. ! Relax for T\n\
1. ! BC salt\n\
1. ! Relax for S\n',
                'b_spillway': '1 0 1 2 ! Bonneville spillway\n\
1. ! Relax for T\n\
0. ! BC salt\n\
1. ! Relax for S\n',
                'b_power2': '0 1 1 2 ! Bonneville power2\n\
1. ! Relax for T\n\
0. ! BC salt\n\
1. ! Relax for S\n'}


OCEAN_BOUNDS = {'pacific': '3 0 0 0 ! Pacific\n',
                'georgia': '3 0 0 0 ! Strait of Georgia\n'}

#-------------------------------------------------------------------------------
# Classes and functions 
#-------------------------------------------------------------------------------
class BCTides(Object):
    " Data and methods to create a bctides.in file "

    def __init__(self, startTime, ndays, grid, ntracers):
        """ 
        Generates a BCTides object to create and write a bctides.in file.
        
        startTime - Run start time as a string in the form of 'mm-dd-yyyy'
        runDays - Length of the run as an integer
        grid - Either the path to the hgrid.gr3 file or specify one 22 or 26
        """
        self.startTime = startTime
        self.loadGrid(grid)
        self.ntracers = ntracers
        self.ndays = ndays


    def loadGrid(self, grid):
        """ Loads a grid into a standard object for parsing """
        # Check if this is a standard grid
        if isinstance(grid,int):
            if grid == 22:
                self.grid = readHGrid("%s/db22/hgrid.gr3" % SCRIPT_PATH)
                self.llPath = "%s/db22/hgrid.ll" % SCRIPT_PATH
            elif grid == 26:
                self.grid = readHGrid("%s/db26/hgrid.gr3" % SCRIPT_PATH)
                self.llPath = "%s/db26/hgrid.ll" % SCRIPT_PATH
            elif grid == 27:
                self.grid = readHGrid("%s/db26c/hgrid.gr3" % SCRIPT_PATH)
                self.llPath = "%s/db26c/hgrid.ll" % SCRIPT_PATH
            else:
                print "I'm unfamiliar with grid %d, please provide a path instead" % grid
        else:
        # Try loading a specified grid
            try:
                self.grid = readHGrid(grid)
            except OSError, e:
                if e.errno == errno.EEXIST:
                    print "The file %s does not exist. Please check path." % grid
                else:
                    print "Problem loading %s, but it appears to exist.  I dunno." % grid
            # Check if hgrid.ll is in the same directory
            if not os.path.isfile(grid[:-3]+'ll'):
                print 'Problem.  I need hgrid.ll to be in the same directory as hgrid.gr3'
                raise
            self.llPath = grid[:-3]+'ll'
                

    def genTidalBoundaries(self):
        """ Generates tidally forced boundaries """
        # Use external programs to generate tidal phase and amplitude on open
        # boundary nodes. Writes to text files.
        genTidalData(self.startTime, self.grid, self.llPath)
        # TODO: Fix this work around. genZ0 needs hgrid.ll and hgrid.gr3 in cwd
        os.system("ln -s %s ./hgrid.ll" % self.llPath)
        hgrid = self.llPath[:-2]+'gr3'
        os.system("ln -s %s ./hgrid.gr3" % hgrid)
        genZ0(self.startTime)

        # Then read that data in from the text files.
        self.readTidalConstituents()

        # Make the nodal factors
        self.genNodalFactors()

        # And then read in the nodal factors
        self.readNodalFactors()

    def genNodalFactors(self):
        """ Generates nodal factor and earth equilibrium for each constituent """

        # Create input file for adcirc tide_fac.f
        f = open('tide_fac.in','w')
        f.write('%d\n' % self.ndays)
        f.write('8 %s' % (time.strftime('%d %m %Y',self.startTime)))
        f.close()

        # Generate the nodal factors by calling tide_fac
        args = ['%s/bin/tide_fac' % SCRIPT_PATH];
        proc = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
        status = proc.wait()    

    def readNodalFactors(self):
        """ 
            Reads in the nodal factors from tide_fac.out 

            self.nfNames - List of names for each nodal factor
            self.nodalFactors - Array of nodal factors and earth equilibrium
                                for each forcing frequency
                self.nodalFactors.shape = (N_CONS_FREQS,2)
        """
        try:
            f = open('tide_fac.out','r')
            nodalFactors = np.zeros((N_CONS_FREQS,2))
            nfNames = []
            # Throw out file header
            for i in xrange(9):         
                tmp = f.readline()
            for i in xrange(N_CONS_FREQS-1):
                tmp = f.readline().rstrip('\n').split()
                nfNames.append(tmp[0])              # Name of constituent
                idx = FREQS_NAMES.index(tmp[0])
                nodalFactors[idx,0] = float(tmp[1])   # Nodal factor
                nodalFactors[idx,1] = float(tmp[2])   # Equilibrium argument
        except:
            raise
        finally:
            f.close()
        # Add Z0 to end of the arrays
        nfNames.append('Z0')
        nodalFactors[0,0] = 0.0
        nodalFactors[0,1] = 0.0
    
        self.nodalFactors = nodalFactors
        self.nfNames = nfNames

    def readTidalConstituents(self):
        """ 
        Reads in tidal constituents at tidally forced nodes from ap_?.dat files.
        Also reads in Z0 from Z0.out which is derived using readssh2b.
        Order corresponds with those in hgrid.gr3.

        self.tidalConstits - List of tidally forced boundaries 
        self.tidalConstits[0] - Array of node phase and amp for boundary 0.
            self.tidalConstits[0].shape = (9,nNodes,2)
        """
        tidalConstits = []
        for ob in glob.glob('ap_?.dat'):
            try:
                f = open(ob,'r')
                f.readline()                                    # Header
                nNodes = int(f.readline().rstrip('\n'))         # Number of nodes on boundary
                f.readline()                                    # Number of constits
                phaseAmp = np.zeros((N_CONS_FREQS,nNodes,2))    # Phase and amplitude at node
                for constit in xrange(N_CONS_FREQS):
                    tmp = float(f.readline().rstrip('\n'))      # Forcing frequency 
                    tmp = f.readline()                          # Constit name
                    for node in xrange(nNodes):
                        tmp = f.readline().rstrip('\n').split()
                        phaseAmp[constit, node, 0] = float(tmp[0])
                        phaseAmp[constit, node, 1] = float(tmp[1])
                tidalConstits.append(phaseAmp)
            except:
                raise
            finally:
                f.close()
        self.tidalConstits = tidalConstits
        self.readZ0()

    def readZ0(self):
        """ 
        Reads in data from Z0.out. MUST BE CALLED AFTER readTidalConstituents()

        self.Z0 - List of arrays of z0 for tidally forced boundaries
        self.Z0[0] - Array of nodal values of z0 along boundary 0
            self.Z0[0].shape = (nNodes)
        """
        Z0 = []
        try:
            f = open("Z0.out", "rb")
            for bound in range(len(self.tidalConstits)):
                nNodes = self.tidalConstits[bound].shape[1]     # Nodes from readTidalConstits
                z0 = np.zeros((nNodes,2))
                f.readline()                                    # Throw out header
                for i in range(nNodes):                         # Get Z0 for all nodes 
                    tmp = f.readline().rstrip('\n').split()
                    z0[i, 0] = float(tmp[0]) 
                    z0[i, 1] = float(tmp[1]) 
                Z0.append(z0)
            f.close()
        except:
            raise

        self.Z0 = Z0

    def writeBoundaries(self, file):
        """ 
        Writes all boundaries to a bctides.in file. Writes all tidally forced
        boundaries first followed by riverine boundaries following precedence.
        It is assumed that the corresponding hgrid.gr3 file follows the same
        standard otherwise there will be problems.

        file - File object for opened bctides.in
        
        """
        file.write('%d ! Number of open boundaries\n' % self.grid.obNodes.nB)
        # Write out ocean boundaries
        idx = 0
        for tb in self.tidalConstits:
            nodes = tb.shape[1]
            idx = idx + 1
            name = eval('self.grid.obNodes.ob_%d_name' % idx)
            file.write('%d %s' % (nodes, OCEAN_BOUNDS[name]))
            nodes = self.tidalConstits[idx-1].shape[1]
            # Write Z0 first because it is derived differently
            file.write('%s\n' % FREQS_NAMES[0])
            for node in xrange(nodes):
                file.write("\t%e\t%e\n" % (self.Z0[idx-1][node][0], 
                                        self.Z0[idx-1][node][1]))
            # Write rest of constituents
            for freq in xrange(1,N_CONS_FREQS):
                file.write('%s\n' % FREQS_NAMES[freq])
                for node in xrange(nodes):
                    file.write("\t%e\t%e\n" % (self.tidalConstits[idx-1][freq][node][0],
                                               self.tidalConstits[idx-1][freq][node][1]))

        # Write out river boundaries
        for ob in range(1, self.grid.obNodes.nB+1):
            var = 'self.grid.obNodes.ob_%d_type' % ob
            if eval(var) is not 'river':
                continue
            name = eval('self.grid.obNodes.ob_%d_name' % ob)
            nodes = eval('self.grid.obNodes.ob_%d_nodes' % ob).shape[0]  
            file.write('%d %s' % (nodes, RIVER_BOUNDS[name]))

    def writeConstFrequencies(self, file):
        """ 
        Writes constituent frequencies to a bctides.in file 

        file - File object of bctides.in
        """
        # Sorted so that they will be in the same order always.
        # Order will differ from what has been done in the past.
        file.write("%d ! Number of tidal constituents\n" % N_CONS_FREQS)
        for fq in xrange(N_CONS_FREQS):
            file.write('%s\n' % FREQS_NAMES[fq])
            file.write('%s %f %f\n' % (FREQS[FREQS_NAMES[fq]], 
                                       self.nodalFactors[fq,0],
                                       self.nodalFactors[fq,1]))

    def writeFile(self):
        """ Writes out gathered data to create a bctides.in file """
        try:
            f = open('bctides.in','w')
            # Header
            str = time.strftime("%m/%d/%Y %H:%M:%S", self.startTime) + " PST"
            f.write(str+'\n')
            f.write("0 40 !ntip and cut off depth earth potential\n")
            # Constituent frequencies, nodal factors, and equilibirum arg
            self.writeConstFrequencies(f)
            self.writeBoundaries(f)
            if self.ntracers > 0:
                self.writeTracers(f)
        except OSError, e:
            if e.errno == errno.EEXIST:
                os.system('mv bctides.in bctides.old')
                print 'Moved old bctides.in to bctides.old'
                self.writeFile()
            else:
                raise
        
        f.close()

    def writeTracers(self, file):
        """ 
            Writes prototype for tracers for bctides.in file 
        
            file - File object of opened bctides.in
        """
        file.write('1 ! itr_met 1:upwind 2:TVD\n')
        file.write('3 ! itrtype(k), k=1, nope_global - 2: timeseries 3: initial value\n')
        for x in xrange(self.grid.obNodes.nB):
            file.write('3\n')
        file.write('0 ! inu_tr ! not valid?')

    def createTemplate(self):
        """ Writes out a generic template for bctides.in """
        # Time to write the file
        try:
            bctides = open('bctides.in', 'wb')
            str = time.strftime("%d/%m%Y %H:%M:%S", self.startTime) + " PST"
            f.write(str)
            bctides.write("ntip - # of contituents in earth tidal potential &\n"
                          "tip_dp - Tidal potential not calculated when depth < tip_dp\n")
            bctides.write("List tidal constituents for each ntip\n")
            bctides.write("nbfr - # of tidal boundary forcing frequencies\n")
            bctides.write("List tidal constituents\n")
            bctides.write("nope - # of open boudary segments\n")
            bctides.write("List elevation, flux, temp, and salinity for each open boundary\n")
            bctides.write("ntracrs - # of tracers\n")
        except:
            print("Error: problem generating bctides.in for test case\n")
            raise
        return

def genBCTides(startTime, ndays, grid, ntracers = 0):
#    if grid == 0:
#        createTemplate()
#        return

    # Create object and get busy child
    bc = BCTides(startTime, ndays, grid, ntracers)

    # Generate boundary conditions
    bc.genTidalBoundaries()

    # Write out the file 
    bc.writeFile()

#-------------------------------------------------------------------------------
# Main 
#-------------------------------------------------------------------------------
if __name__ == '__main__':
    usage = ("%prog [startTime mm-dd-yyyy]) [rundays] and specify path to hgrid.gr3 "
             "via -p or specify a grid (22,26) via -g")

    parser = OptionParser(usage=usage)
    parser.add_option("-p", "--path", action="store", type="string",
                      dest="path", help="Path to hgrid.gr3")
    parser.add_option("-g", "--grid", action="store", type="int",
                      dest="grid", help="Grid number [22, 26]")
    parser.add_option("-t", "--ntracers", action="store", type="int",
                      dest="ntracers", help="Number of tracers for run")

    parser.set_defaults(path='')
    parser.set_defaults(grid=999)
    parser.set_defaults(ntracers=0)

    (options, args) = parser.parse_args()
    grid = options.grid
    path = options.path
    ntracers = options.ntracers
    if (grid == 999 and path == '') or (grid != 999 and path != ''):
        parser.error("Apply either a path to a grid or a grid number")
    if path != '':
        grid = path

# grab start and number of days 
    if len(args) != 2:
        parser.error("Incorrect number of arguments")
    else:
        start = args[0]
        ndays = int(args[1])

# deal with start days
    (month, day, year) = start.split('-')
    startDay   = int(day)
    startMonth = int(month)
    startYear  = int(year) 
    startTime = datetime.datetime(startYear, startMonth, startDay)

    genBCTides(startTime, ndays, grid, ntracers)
