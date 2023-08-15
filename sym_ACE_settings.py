# global settings for the sym_ACE library
import sys
import pathlib

# store CG and wigner files in top_dir/lib/

import gen_labels
topfile = gen_labels.__file__
top_dir = topfile.split('/ge')[0]
lib_path = '%s/lib' % top_dir

# store CG and wigner files in current directory:

#lib_path = str(pathlib.Path(__name__).parent.resolve())

#-------------------------------
# Coupling coefficient settings
#-------------------------------
ranks = range(1,6)
# lmax per rank
lmax_dict_G = {1:0 , 2:8, 3:6, 4:4, 5:2}
# lmax for underlying traditional wigner-3j symbols (will impose limits on lmax_dict_G)
lmax_traditional=12
# Flag to generate library of CG coefficients as well
cglib = True
# Flag to generate generalized wigner coefficients for L_R = 1 reduced reps (vector-like descriptors)
gen_LR1 = False
