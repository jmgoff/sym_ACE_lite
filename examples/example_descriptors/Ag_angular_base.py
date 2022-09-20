from descriptor_calc_local.RI_base import *

#-----------------------------------------------------------
# INPUT and definition of the model
#-----------------------------------------------------------
ranks = range(1 , 6 + 1)                                    # ranks of basis functions to be evaluated
elements = ['Ag']                                 # list of elements
rc = 20.5                                              # radial cutoff (angstroms) for radial basis
lmbda=5.0                                             # exponential factor for sampling of r << r_c in g_k(r)
nradbase=16                                           # maximum k in g_k expansion of R_nl (expansion representation of radial basis not implemented yet) this currently just specifies the number of rank 1 (2-body) functions
lmax_dict = {1:0,2:2,3:2,4:2,5:2,6:2}                 # maximum l quantum number in basis functions per rank: {rank:lmax(rank)} (note lmax of rank 1 has no effect)
nradmax_dict = {1:nradbase, 2: 2, 3:2, 4:2, 5:2, 6:2}   # maximum n quantum number in basis functions per rank: {rank:nmax(rank)}
i_d = 0                                               # index supplied to identify the atomic configuration (useful when iterating over many configurations)


# flatten the descriptor list for function input
nus = ['1,1,1,1']

# example structure
from ase import Atoms,Atom
atoms = Atoms(['Ag']*10)
np.random.seed(55689)
positions = np.random.rand(10,3)
positions = 5.0 * positions
atoms.set_positions(positions)
atoms.set_pbc(False)


#-----------------------------------------------------------
# RUNNING the analytical ACE code to get per-atom descriptor
#-----------------------------------------------------------

# input for descriptor calculator is expected in dictionary format (in order to use JSON input files later)
arg = {'id':i_d,'atoms':atoms,'elements':elements,'rc':rc,'nradbase':nradbase,'nradmax':nradmax_dict,'lmax':lmax_dict,'lmbda':lmbda,'rank':ranks,'nus':nus, 'L_R':0, 'M_R':0}

B_per_atom = atoms_RI(arg)

import json
with open('symmetric_ag.json','w') as writejson:
	json.dump(B_per_atom,writejson, sort_keys=False, indent=2)
