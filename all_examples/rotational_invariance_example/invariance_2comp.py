from equivariant_calc.chem_B import *
from gen_labels import *
from scipy.spatial.transform import Rotation as R

#-----------------------------------------------------------
# INPUT
#-----------------------------------------------------------
ranks = range(1,5)                    # ranks of basis functions to be evaluated
rc = 22.5                          # radial cutoff (angstroms) for radial basis
lmbda=5.0                         # exponential factor for sampling of r << r_c in g_k(r)
nradbase=16                       # maximum k in g_k expansion of R_nl (expansion representation of radial basis not implemented yet)
lmax_dict = {1:6,2:2,3:2,4:2}         # maximum l quantum number in basis functions per rank: {rank:lmax(rank)} 
nradmax_dict = {1:nradbase, 2: 2, 3:1, 4:1}   # maximum n quantum number in basis functions per rank: {rank:nmax(rank)}
i_d = 0                           # index supplied to identify the atomic configuration (useful when iterating over many configurations)


from ase.io import read,write
#atoms = read('Cu_4.xyz')                  # ASE atoms object
#atoms = read('CuAg.cif')                  # ASE atoms object
atoms = read('random_CuAg.xyz')                  # ASE atoms object



#example for generating valid l_1,l_2, ... l_N labels for N spherical harmonics to be reduced
lrng = [1]
n_to_reduce = 4
L_R = 0
M_R= 0

# example for generating couplings used to reduce spherical harmonics
from clebsch_couple import *
ldict = {n_to_reduce:max(lrng)}
try:
    #with open('cg_LR_0_r4_lmax3.pickle','rb') as handle:
    #with open('cg_LR_12_r4_lmax3.pickle','rb') as handle:
    with open('cg_LR_%d_r4_lmax%d.pickle' %(L_R,max(lrng)),'rb') as handle:
    #with open('cg_LR_2_r4_lmax3.pickle','rb') as handle:
    #with open('cg_LR_4_r4_lmax3.pickle','rb') as handle:
    #with open('cg_LR_6_r4_lmax3.pickle','rb') as handle:
    #with open('cg_LR_8_r4_lmax3.pickle','rb') as handle:
        ccs = pickle.load(handle)
except FileNotFoundError:
    ccs = get_cg_coupling(ldict,L_R=L_R)
    print (ccs)

#-----------------------------------------------------------
# RUNNING the analytical ACE code to get per-atom descriptor
# values in an atomic configuration (can easily be summed
# to obtain dE/dc entries in the A matrix)
#-----------------------------------------------------------

# print the number of descriptors being calculated

# flatten the list for function input
#nus = [item for sublist in ranked_nus for item in sublist]
nus = ['0_0,0,0,0,1,1,1,1,1,1,1,1_2-2']
mus = {nu:[] for nu in nus}
#chems = ['Cu','Ag']
chems = ['Ag','Cu']
for nu in nus:
    vcs = ind_vec(chems,4)
    for mustr in vcs:
        mus[nu].append(mustr)
elements = ['Ag','Cu']
descriptor_type = 'chemical_radial_angular'
# input for descriptor calculator is expected in dictionary format (in order to use JSON input files later)
#arg = {'id':i_d,'atoms':atoms,'rc':rc,'nradbase':nradbase,'nradmax':nradmax_dict,'lmax':lmax_dict,'lmbda':lmbda,'rank':ranks,'nus':nus,'mus':mus}
arg = {'id':i_d,'atoms':atoms,'elements':elements,'rc':rc,'nradbase':nradbase,'nradmax':nradmax_dict,'lmax':lmax_dict,'lmbda':lmbda,'rank':ranks,'nus':nus,'ccs':ccs, 'descriptor_type':descriptor_type, 'L_R':0, 'M_R':0}

B_per_atom = atoms_eB(arg)
B_per_atom2 = atoms_eB(arg, **{'transformation':R.from_rotvec(np.pi/2 * np.random.rand(3))})
print (B_per_atom,B_per_atom2)
import json
with open('random_B.json','w') as writejson:
    json.dump(B_per_atom,writejson, sort_keys=True, indent=2)
with open('random_B_rot.json','w') as writejson:
    json.dump(B_per_atom2,writejson, sort_keys=True, indent=2)
print ('atom_index, nl, mu, descriptor, descriptor_rotated, difference')
for atind in range(len(atoms)):
    for nu in nus:
        for mu in mus[nu]:
            try:
                v1 = B_per_atom[atind][nu]
            except KeyError:
                v1 = 0.0
            try:
                v2 = B_per_atom2[atind][nu]
            except KeyError:
                v2 = 0.0
            diff = v2-v1
            print (atind,nu,mu,v1,v2, diff)
