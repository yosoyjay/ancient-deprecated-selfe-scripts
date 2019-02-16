import sys
sys.path.append('/home/hyde/ourestuaries/pythonlib/cmop')
from meccs_file_util import *

input_file = "/home/users/hyde/transects.txt"
out_base = "/home/hyde/test_tran"

cs_to_csvs(input_file, out_base)

#input_file = "/home/hyde/public_html/ourestuaries/shapes_data/hyde/tilc/cross_sections_test.txt"
#test_file = "/home/hyde/public_html/ourestuaries/shapes_data/hyde/tilc/cross_sections_test_out.txt"
#old_file = "/home/hyde/public_html/ourestuaries/shapes_data/hyde/tilc/cross_sections.txt"
#new_old_file = "/home/hyde/public_html/ourestuaries/shapes_data/hyde/tilc/cross_sections_old_out.txt"
#
#(pts, names, rivers) = read_cs_file(input_file, 1)
#
#write_cs_file(test_file, pts, names, rivers)
#
#(pts2, names2, rivers2) = read_cs_file(test_file, 1)
#
#(pts3, names3, rivers3) = read_cs_file(old_file, 1)
#
#write_cs_file(new_old_file, pts3, names3, rivers3)
#
#(pts4, names4, rivers4) = read_cs_file(new_old_file, 1)
#
#pdb.set_trace()
#
#dbg = 1
