# global settings for the sym_ACE library

# Name of the directory the repo is stored in:
pkg_name = 'sym_ACE_lite'

#-------------------------------
# Coupling coefficient settings
#-------------------------------
ranks = range(1,7)
# lmax per rank
lmax_dict_G = {1:0 , 2:8, 3:6, 4:4, 5:3, 6:2}
# lmax for underlying traditional wigner-3j symbols (will impose limits on lmax_dict_G)
lmax_traditional=10
# Flag to generate library of CG coefficients as well
cglib = False
# Flag to generate generalized wigner coefficients for L_R = 1 reduced reps (vector-like descriptors)
gen_LR1 = False
