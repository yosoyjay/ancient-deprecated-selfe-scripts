#!/usr/bin/python
# This script generates hgrid.new and depth_diff.gr3 following
# Joe Cho's guide.
#
# jlopez 3/12/2010
#

# Some constants like directories and such
homePath  = "/home/workspace/project/lopezj"
binPath  = "/home/workspace/project/lopezj/src/selfe_run_prep/src"

import os
import sys
sys.path.append("%s/bin/pexpect/" % homePath)
import pexpect

# Run deep_smooth program
os.system("%s/src/deep_smooth")
