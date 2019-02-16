#!/usr/local/bin/python
"""This script sets up the folders and files for a SELFE run. 

The run is set up assuming "regularly" used values for a hindcast.  Please see
help for more details about changes to be made via CLI arguments.

Note: Some modules are imported immediately before the functions are called
due to dependencies fulfilled only at that point and poor planning on my
part.
"""
import os 
import sys
import shutil 
import datetime 
import subprocess 
from optparse import OptionParser
from change_gr3 import change_gr3

# constant stuff and imported functions
SCRIPT_PATH = '/home/workspace/users/lopezj/scripts/selfe_setup/'
sys.path.append(SCRIPT_PATH)

################################################################################
#
# __createDirs - Creates directories for a run 
#
################################################################################

def createDirs(startDate, grid, run_name):
    """ Creates run directories
    """

    try:
        if run_name == '':
            run_dir = "%i-%s-%s_db%i_%s_days" % (startDate[0],                   
                str(startDate[1]).zfill(2), str(startDate[2]).zfill(2), grid,   
                str(n_days).zfill(2))
        else:
            run_dir = run_name
        os.mkdir("%s" % run_dir) 
        os.mkdir("%s/run" % run_dir) 
        os.mkdir("%s/run/outputs" % run_dir) 
        os.mkdir("%s/post" % run_dir) 
        os.mkdir("%s/log" % run_dir)
        os.mkdir("%s/setup" % run_dir)
        os.system("mv setup_params.log %s/log " % run_dir)
        os.chdir("%s/run" % run_dir) 
    except OSError:
        print("Error: problem creating run directories")
        print("Make sure a directory name %s does not already exist" % run_dir)
        raise

################################################################################
#
# __createLinks- Creates links for pelfe, autocombine.pl, and selfe_run 
#
################################################################################

def createLinks(startDate, grid, email):
    try:
        os.system('cp -s %s/../../bin/pelfe ./pelfe' % SCRIPT_PATH)
#       shutil.copy("%s/../../bin/pelfe" % SCRIPT_PATH, "./pelfe")
    except:
        print("Error: Problem linking pelfe to run dir")
        raise

    try:
        shutil.copy('%s/selfe_run' % SCRIPT_PATH, './selfe_run')
    except:
        print("Error: Problem copying selfe_run to run dir")
        raise

    try:
        shutil.copy("%s/autocombine.pl" % SCRIPT_PATH, "./autocombine.pl")
    except:
        print("Error: Problem copying autocombine.pl to run dir")
        raise

# If user wants an email message about run, do it.
# Maybe create an object for the selfe_run script 
    if email:
        try:
            _selfe_run = open("selfe_run~", "wb")
            _cmd_args = ["sed", "/#$ -S/{p;s/.*/#$ -M %s/;p;s/.*/#$ -m ae/;}"
                          % email, "selfe_run"]
            subprocess.check_call(_cmd_args, stdout = _selfe_run)
            _selfe_run.close()
            _cmd_args = ["mv", "-f", "selfe_run~", "selfe_run"]
            subprocess.check_call(_cmd_args)
        except:
            print("Error: problem editing selfe_run.  Check file for errors.")
            raise
# Name the run based on db and date
    try:
        _selfe_run = open("selfe_run~", "wb")
        _cmd_args = ["sed", "/#$ -S/{p;s/.*/#$ -N db%i_%i%i/;}" % (grid, 
                    startDate[1], startDate[0]), "selfe_run"]
        subprocess.check_call(_cmd_args, stdout = _selfe_run)
        _selfe_run.close()
        _cmd_args = ["mv", "-f", "selfe_run~", "selfe_run"]
        subprocess.check_call(_cmd_args)
    except:
        print("Error: problem editing selfe_run.  Check file for errors.")
        raise

################################################################################
#
# Generate a gr3 file 
#
################################################################################

def __genGR3(file_name):
    print "Generating %s\n" % file_name
    gr3_files = {"diffmin.gr3":0.000001,
                 "diffmax.gr3":0.01,
                 "xlsc.gr3":0.5,
                 "drag.gr3":0.0045,
                 "interpol.gr3":2}
    value = raw_input("Typical value %f\n"
                      "Enter a value for %s:" % (gr3_files[file_name], file_name))

    change_gr3('hgrid.gr3', file_name, -999, float(value), float(value))



################################################################################
#
# Creates a very generic test step up directory structure with essential files 
#
################################################################################

def testSetup(startDate, n_days, grid, grid_path, n_tracers, email='', time_step=90,
              run_name=''):

    endDate   = startDate + datetime.timedelta(days=(n_days-1))
    startDate = startDate.timetuple()       
    endDate   = endDate.timetuple()
    # If year is > 2005 use NAM as sflux source otherwise use NARR
    if startDate[0] >= 2005:
        sflux_src = 1
    else:
        sflux_src = 3


################################################################################
#
# Make directory structure
#
################################################################################

    print 'Creating run directories'
    createDirs(startDate, grid, run_name)

################################################################################
#
# Check hgrid.gr3 and generate *.gr3 files 
#
################################################################################
    if os.path.exists(grid_path):
        shutil.copy(grid_path, "./hgrid.gr3")
    else:
        sys.exit("Error: Unable to find %s\n"
                 "Please check the path.\n")

    print 'Generating gr3 files for test run'
    __genGR3('diffmax.gr3')
    __genGR3('diffmin.gr3')
    __genGR3('interpol.gr3')
    __genGR3('xlsc.gr3')
    __genGR3('drag.gr3')

################################################################################
#
# Generate bctides.in 
#
################################################################################

    print '\nGenerating bctides.in\n'
    from genBCTides import genBCTides
    genBCTides(startDate, n_days, grid_path, n_tracers)

################################################################################
#
# Copy param.in and change number of days for run and turn on writing 
# tracer output files.  Only turns on two.
#
################################################################################
    try:
        shutil.copy('%s/param.in' % SCRIPT_PATH, './param.in')
    except:
        print("Error: problem copying param.in from %s to run dir" % SCRIPT_PATH)
        raise

    try:
        _cmd_args = ["sed", "-i", "151s/ X / %i /" % n_days, "param.in"]
        subprocess.check_call(_cmd_args)
# if the run time step is not 90, then changes must be made to ihfskip, 
# hotout_write, and nspool to make sure they are produced at the same rate.
# ihfskip one a day and nspool every 15 minutes.
# ihfskip and hotout_write typically have the same value.
        if time_step != 90:
            ihfskip = 86400/time_step
            nspool = 900/time_step
            _cmd_args = ["sed", "-i", "262s/960/%i /" % ihfskip, "param.in"]
            subprocess.check_call(_cmd_args)
            _cmd_args = ["sed", "-i", "301s/960/%i /" % ihfskip, "param.in"]
            subprocess.check_call(_cmd_args)
            _cmd_args = ["sed", "-i", "261s/10/%i /" % nspool, "param.in"]
            subprocess.check_call(_cmd_args)
            _cmd_args = ["sed", "-i", "154s/90./%i /" % time_step, "param.in"]
            subprocess.check_call(_cmd_args)
        if n_tracers > 0:
            _cmd_args = ["sed", "-i", "24s/ 0 / %i /" % n_tracers, "param.in"]
            subprocess.check_call(_cmd_args)
# TODO: Okay, this is stupid. Make a param.in parser. 
            for i in range(n_tracers):
                _line_number = 290
                _cmd_args = ["sed", "-i", "%is/ 0/ 1/" % _line_number, 
                             "param.in"]
                subprocess.check_call(_cmd_args)
                _line_number += 1
                if i==2:
                    print("You must create tracer ouput files in param.in if"
                          " you are using more than 2 tracers.\n")
                    break
    except:
        print("Error: problem editing param.in. Check file for errors.")
        raise

################################################################################
#
# Create links 
#
################################################################################

    print "Copy or linking to pelfe, autocombine.pl, and qsub script (selfe_run)"
    createLinks(startDate, grid, email) 

    print ("Finished generated run directories for a test run.\n"
           "Please check param.in and bctides.in for appropriate values")





def selfeSetup(startDate, n_days, grid, n_tracers, email='', time_step=90, 
               hot_start=1, ncom_version=0, run_name='', th_exist=0):
    """ Main function that essentially automates Joe Cho's guide to running 
        Selfe. 
    """

    endDate   = startDate + datetime.timedelta(days=(n_days-1))
    startDate = startDate.timetuple()       
    endDate   = endDate.timetuple()
    # Use NARR 3, NAM is 1
    sflux_src = 1

################################################################################
#
# Make directory structure
#
################################################################################

    print 'Creating run directories'
    createDirs(startDate, grid, run_name)

################################################################################
#
# Link *.gr3 files if using standard grid or generate those files if using a
# test grid.
#
# Grid 15 is hi-res.
# Grid 27 is grid 26 cut at Beaver Army.
# Grid else is hi-res 26. 
#
################################################################################

    print '\nCopying gr3 files for db%i' % grid
    if(grid == 15):
        ret_code = os.system('cp %s/hi_res/* ./' % SCRIPT_PATH)
        if ret_code != 0:
            print("Error: problem creating symlinks of grid files to run directory")
            raise IOError
    elif(grid == 22):
        ret_code = os.system('cp %s/db22/* ./' % SCRIPT_PATH)
        if ret_code != 0:
            print("Error: problem creating symlinks of grid files to run directory")
            raise IOError
    elif(grid == 26):
        ret_code = os.system('cp %s/db26/* ./' % SCRIPT_PATH)
        if ret_code != 0:
            print("Error: problem creating symlinks of grid files to run directory")
            raise IOError
    elif(grid == 27):
        ret_code = os.system('cp %s/db26c/* ./' % SCRIPT_PATH)
        if ret_code != 0:
            print("Error: problem creating symlinks of grid files to run directory")
            raise IOError
    else: 
# High resolution version of 26, set up grid stuff then the rest is based on 26.
        #os.symlink('%s/db22/*' % SCRIPT_PATH, './')
        ret_code = os.system('cp %s/db26_div4/* ./' % SCRIPT_PATH)
        if ret_code != 0:
            print("Error: problem creating symlinks of grid files to run directory")
            raise IOError
        grid = 26

################################################################################
#
# Generate sflux - using nam (1) as a source 
#
################################################################################

    print "\nGenerating sflux\n"
    command = "%s/make_sflux_links.csh" % SCRIPT_PATH
    try:
        sflux = open("../log/sflux.log","w")
        return_code = subprocess.call([command, str(sflux_src), 
            str(startDate[0]), str(startDate[1]), str(startDate[2]),
            str(endDate[0]), str(endDate[1]), str(endDate[2])], stdout=sflux)
        if return_code != 0:
            sys.exit("Error: Problem generating sflux.")
    except:
        print("Error: Problem generating sflux.")
        raise

################################################################################
#
# Generate tidal data
#
################################################################################

    print '\nGenerating tidal data\n'
    from genTidalData import genTidalData
    genTidalData(startDate, grid)

################################################################################
#
# Generate Z0 (sea-surface) values 
#
################################################################################
    print '\nGenerating Z0 values\n'
    from genZ0 import genZ0
    genZ0(startDate, grid)

################################################################################
#
# Generate bctides.in 
#
################################################################################

    print '\nGenerating bctides.in\n'
    from genBCTides   import genBCTides
    genBCTides(startDate, grid)

################################################################################
#
# Copy param.in and change number of days for run and turn on writing 
# tracer output files.  Only turns on two.
#
################################################################################
    try:
        shutil.copy('%s/param.in' % SCRIPT_PATH, './param.in')
    except:
        print("Error: problem copying param.in from %s to run dir" % SCRIPT_PATH)
        raise

    try:
        _cmd_args = ["sed", "-i", "151s/ X / %i /" % n_days, "param.in"]
        subprocess.check_call(_cmd_args)
# if the run time step is not 90, then changes must be made to ihfskip, 
# hotout_write, and nspool to make sure they are produced at the same rate.
# ihfskip one a day and nspool every 15 minutes.
# ihfskip and hotout_write typically have the same value.
        if time_step != 90:
            ihfskip = 86400/time_step
            nspool = 900/time_step
            _cmd_args = ["sed", "-i", "262s/960/%i /" % ihfskip, "param.in"]
            subprocess.check_call(_cmd_args)
            _cmd_args = ["sed", "-i", "301s/960/%i /" % ihfskip, "param.in"]
            subprocess.check_call(_cmd_args)
            _cmd_args = ["sed", "-i", "261s/10/%i /" % nspool, "param.in"]
            subprocess.check_call(_cmd_args)
            _cmd_args = ["sed", "-i", "154s/90./%i /" % time_step, "param.in"]
            subprocess.check_call(_cmd_args)
        if n_tracers > 0:
            _cmd_args = ["sed", "-i", "24s/ 0 / %i /" % n_tracers, "param.in"]
            subprocess.check_call(_cmd_args)
# TODO: Okay, this is stupid.  Fix this. 
            for i in range(n_tracers):
                _line_number = 290
                _cmd_args = ["sed", "-i", "%is/ 0/ 1/" % _line_number, 
                             "param.in"]
                subprocess.check_call(_cmd_args)
                _line_number += 1
                if i==2:
                    print("You must create tracer output files in param.in if"
                          " you are using more than 2 tracers.\n")
                    break
    except:
        print("Error: problem editing param.in. Check file for errors.")
        raise

################################################################################
# 
# Create links to pelfe binary, qsub script (selfe_run), and MPI combine script
#
################################################################################

    try:
#       os.system('ln -s %s/../../bin/pelfe ./pelfe' % SCRIPT_PATH)
        shutil.copy("%s/../../bin/pelfe" % SCRIPT_PATH, "./pelfe")
    except:
        print("Error: Problem linking pelfe to run dir")
        raise

    try:
        shutil.copy('%s/selfe_run' % SCRIPT_PATH, './selfe_run')
    except:
        print("Error: Problem copying selfe_run to run dir")
        raise

    try:
        shutil.copy("%s/autocombine.pl" % SCRIPT_PATH, "./autocombine.pl")
    except:
        print("Error: Problem copying autocombine.pl to run dir")
        raise

# If user wants an email message about run, do it.
# Maybe create an object for the selfe_run script 
    if email:
        try:
            _selfe_run = open("selfe_run~", "wb")
            _cmd_args = ["sed", "/#$ -S/{p;s/.*/#$ -M %s/;p;s/.*/#$ -m ae/;}"
                          % email, "selfe_run"]
            subprocess.check_call(_cmd_args, stdout = _selfe_run)
            _selfe_run.close()
            _cmd_args = ["mv", "-f", "selfe_run~", "selfe_run"]
            subprocess.check_call(_cmd_args)
        except:
            print("Error: problem editing selfe_run.  Check file for errors.")
            raise
# Name the run based on db and date
    try:
        _selfe_run = open("selfe_run~", "wb")
        _cmd_args = ["sed", "/#$ -S/{p;s/.*/#$ -N db%i_%i%i/;}" % (grid, 
                    startMonth, startYear), "selfe_run"]
        subprocess.check_call(_cmd_args, stdout = _selfe_run)
        _selfe_run.close()
        _cmd_args = ["mv", "-f", "selfe_run~", "selfe_run"]
        subprocess.check_call(_cmd_args)
    except:
        print("Error: problem editing selfe_run.  Check file for errors.")
        raise


################################################################################
#
# Generate temp.th and flux.th first becase temp.th is needed for hotstart.in
#
################################################################################

    print '\nGenertating flux.th and temp.th\n'
    from genFluxTemp  import genFluxTemp
    genFluxTemp(startDate, n_days, grid, time_step, th_exist)

################################################################################
#
# Generate hotstart.in and nudging files 
#
################################################################################

    print hot_start
    if hot_start == 1:
        from genHotstart import genHotstart
        genHotstart(startDate, ncom_version, grid, n_tracers)
    else:
# Create date.in file that will be used for nudging files
        try:
            datein = open("date.in", "w")
            try:
                if(grid==15):
                    datein.write('1 1 1\n')     # 1 open boundary of type 4
                else:
                    datein.write('1 1 2\n')     # 2 open boundary of type 4?
                datein.write('%i %i %i 0\n' % 
                    (startDate[0], startDate[1], startDate[2]))    
                datein.write('0\n')             # ireduce
                datein.write('Selfe\n')         # selfe or elcirc
                datein.write('1.e6 33.\n')      # ht(depth), salinitymin 
                datein.write('%s\n' % n_tracers)# number of tracers
            finally:
                datein.close()
        except IOError:
            print("error: unable to open date.in for writing\n")
            raise

# Change date.in and call readncom to generate nudging files from NCOM
    print '\nGenerating nudging files\n'
    try:
        _cmd_args = ["sed", "-i", "2s/.$/%i/" % n_days, "date.in"]  
        subprocess.check_call(_cmd_args)
    except:
        print("Error: Problem editing date.in for nudging files")
        raise
    try:
        _cmd_args = "%s/bin/readncom8b" % SCRIPT_PATH
        subprocess.check_call(_cmd_args)
    except:
        print("Error: problem generating nudging files")
        raise


################################################################################
#
# Save set up files in directory for archiving and to facilitate troubleshooting 
#
################################################################################

    files = ['ap_1.dat', 'ap_2.dat', 'boundary.new', 'date.in', 'intel_deg', 
            'nos8', 'nos8_0.new', 'nos8_1.new', 'ssh.in', 'table.new',
            'tide_colu.com', 'tide_colu.out', 'Z0.out']

    for file in files:
        shutil.move(file, "../setup/")

    print ("Finished. Baseline run setup, specific implementation details not" 
           " handled by the script and must be made manually.\n")

################################################################################
#
# Main - Check CLI arguments 
#
################################################################################

if __name__ == '__main__':
# Create a parser for command line args and options
    usage = ("Usage: %prog [begin_date mm-dd-yyyy] [number of days in run] "
             "[grid [15(hi-res), 22, 26, 27(26c), 27(26 hi-res)]] [options]")

    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--tracer", action="store", type="int",
                      dest="n_tracers", help="number of tracers")
    parser.add_option("-m", "--mail", action="store", type="string",
                      dest="email_addr", help="email to send run message")
    parser.add_option("-d", "--delta", action="store", type="int",
                      dest="time_step", help="time step of run")
    parser.add_option("-H", "--hotstart", action="store", type="int",
                      dest="hot_start", help="gen hotstart.in")
    parser.add_option("-N", "--ncom", action="store", type="int",
                      dest="ncom", help="NCOM ic version")
    parser.add_option("-g", "--grid", action="store", type="string",
                      dest="test_grid", help="Test grid path")
    parser.add_option("-n", "--name", action="store", type="string",
                      dest="run_name", help="Name of the run")
    parser.add_option("-T", "--th", action="store", type="string",
                      dest="th_exist", help="temp.th and flux.th exist")

    parser.set_defaults(n_tracers=0)
    parser.set_defaults(time_step=90)
    parser.set_defaults(hot_start=1)
    parser.set_defaults(ncom=0)
    parser.set_defaults(test_grid='')
    parser.set_defaults(run_name='')
    parser.set_defaults(th_exist=0)

    (options, args) = parser.parse_args()
    n_tracers = options.n_tracers
    time_step = options.time_step
    hot_start = options.hot_start
    grid_path = options.test_grid
    ncom = options.ncom
    run_name = options.run_name 
    th_exist = options.th_exist

# If email is not used, send it as empty string.  
    if options.email_addr:
        email = options.email_addr
    else:
        email = ''

# If a test grid is specified, set grid to a magic number 0 and check path
    if grid_path: 
        grid = 0 
        if os.path.exists(grid_path):
            grid_path = os.path.abspath(grid_path)
        else:
            sys.exit("Error: Unable to locate %s" % grid_path)

# grab start date, number of days, and db(grid) to use
    if len(args) != 3:
        parser.error("Incorrect number of arguments")
    else:
        start  = args[0]
        n_days = int(args[1])
        grid   = int(args[2])

# if not a known grid die. 
    grids = [15, 22, 26, 27, 28, 0]
    if grid not in grids:
        sys.exit("Grid options are: 15(hi-res), 22, 26, 27(26c), "
                 "28(hi-res 26), or 0 to indicate user supplied grid")

# deal with start and ending days
    (month, day, year) = start.split('-')
    startDay   = int(day)
    startMonth = int(month)
    startYear  = int(year) 

# python style dates
    startDate = datetime.date(startYear, startMonth, startDay)

################################################################################
#
# Write the parameters to a log file to ensure there is a record of the run 
# parameters.
#
################################################################################
    try:
        logFile = open('setup_params.log', 'w')
        logFile.write('Start Date: %d-%d-%d\n' % (startMonth, startDay, startYear))
        logFile.write('Days in run: %d\n' % n_days)
        logFile.write('Grid used: %d\n' % grid)
        if grid==0:
            logFile.write('Test grid path: %s\n' % grid_path)
        logFile.write('Number of tracers: %d\n' % n_tracers)
        logFile.write('Time step: %d\n' % time_step)
        logFile.write('Hotstart created?: %d\n' % hot_start)
        logFile.write('NCOM Scaling factor: %d\n' % ncom)
        logFile.write('Email sent to: %s\n' % email)
    except:
        print("Problem creating setup_params.log file")
        raise
    finally:
        logFile.close()

    if grid!=0:
        selfeSetup(startDate, n_days, grid, n_tracers, email, time_step, 
            hot_start, ncom, run_name, th_exist)
    else:
        testSetup(startDate, n_days, grid, grid_path, n_tracers, email, time_step, run_name)
