import sys
sys.path.append('/home/hyde/ourestuaries/pythonlib/cmop')
from numpy import *
from meccs_file_util import *

if len(sys.argv) < 3:
    sys.exit("usage: python gen_csv.py [cs_file] [output_base]\n")

cs_file = sys.argv[1]
base_dir = sys.argv[2]

cs_to_csvs(cs_file, base_dir)

