import sys
sys.path.append('/home/hyde/ourestuaries/pythonlib/cmop')
from numpy import *
from gr import *
from ob import *
import shape_procs
from model_util import *
#import pdb
from model_reader import mreader
import meccs_file_util
#import sys
from pylab import save, load
import os.path
import pdb

fname = "/home/hyde/public_html/ourestuaries/shapes_data/hyde/yaqalsea/transects.ob"

tob = Ob()
tob.load_ob(fname)

fname = "/home/hyde/public_html/ourestuaries/shapes_data/hyde/yaqalsea/cross_sections.ob"

cob = Ob()
cob.load_ob(fname)

pdb.set_trace()

dbg = 1
