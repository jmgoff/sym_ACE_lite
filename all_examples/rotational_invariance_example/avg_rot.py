#from wigner_couple import *
from clebsch_couple import *
#from get_sym import *
from equivariant_calc.chem_B import *
from scipy.spatial.transform import Rotation as R

from ase.build import bulk

import json

from ase import Atoms,Atom
from ase.io import read,write

import numpy as np
#-----------------------------------------------------------
# INPUT
#-----------------------------------------------------------
ranks = [1,4]				    # ranks of basis functions to be evaluated
rc = 20.					      # radial cutoff (angstroms) for radial basis
lmbda=5.0					     # exponential factor for sampling of r << r_c in g_k(r)
#lmbda=0.005					     # exponential factor for sampling of r << r_c in g_k(r)
nradbase=16					   # maximum k in g_k expansion of R_nl (expansion representation of radial basis not implemented yet)
lmax_dict =    {1:0,        2: 2, 3:2, 4:3,5:4,6:1}		 # maximum l quantum number in basis functions per rank: {rank:lmax(rank)}
nradmax_dict = {1:nradbase, 2: 2, 3:2, 4:5,5:3,6:1}   # maximum n quantum number in basis functions per rank: {rank:nmax(rank)}

elements = ['Mo'] #should be in 

import sys

i_d = int(sys.argv[1])


#np.random.seed(123)
np.random.seed(i_d)
r = R.from_rotvec(np.pi/2 * np.random.rand(3))


unsc = np.random.uniform(-1,1,(20,3))
pos = unsc*7

pos_rot = np.matmul(r.as_matrix(),pos.T)
pos_rot = pos_rot.T

pos_rot = list(pos_rot)
pos_rot = [np.array([0.,0.,0.])] + pos_rot

atoms = Atoms(symbols=['Mo']*len(pos_rot))
atoms.positions = pos_rot
atoms.set_cell([[14,0,0],[0,14,0],[0,0,14]])
write('Mo_rot.data',atoms,format='lammps-data')


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
	#store them for later so they don't need to be recalculated
	store_generalized(ccs, coupling_type='cg',L_R=L_R)



orb_nls = [ #'cg_LR_0_r4_lmax3.pickle'
'0_0,0,0,0,1,1,1,1,1,1,1,1_0-0' ,
'0_0,0,0,0,1,1,1,1,1,1,1,1_2-2' ,
]

descriptor_type = 'chemical_radial_angular'
for orb_nl in orb_nls:
	nus = [orb_nl]
	#mu0munlL,rank = get_nu_vectors(orb_nl)
	mu0,mu,n,l,inter = get_mu_n_l(orb_nl,return_L=True)
	rank = get_mu_nu_rank(orb_nl)
	nodes,remainder = tree(l)
	base_node_inters = {node:get_intermediates_w(l[node[0]],l[node[1]]) for node in nodes}
	mstrs = get_ms(l)

	arg = {'id':i_d,'atoms':atoms,'elements':elements,'rc':rc,'nradbase':nradbase,'nradmax':nradmax_dict,'lmax':lmax_dict,'lmbda':lmbda,'rank':ranks,'nus':nus, 'descriptor_type':descriptor_type, 'L_R':0, 'M_R':0}
	lint = inter
	lint_str_list = ['%d'] * len(lint)
	lint_str = '-'.join(b for b in lint_str_list) % tuple(lint)
	nstr_list = ['%d']*len(n)
	nstr = ''.join(b for b in nstr_list) % tuple(n)
	lstr_list = ['%d']*len(l)
	lstr = ''.join(b for b in lstr_list) % tuple(l)
	lstrp = ','.join(b for b in lstr_list) % tuple(l)
	print(lint)
	if not os.path.isfile('struct%s_%s_%d_%s_rot.json' % (nstr,lstr,i_d,lint_str)):
		#B_per_atom = atoms_eB(arg,**{'ccs':ccs[M_R][lstrp][inter]})
		B_per_atom = atoms_eB(arg,**{'ccs':ccs[M_R]})
		ats = list(B_per_atom.keys())
		print (B_per_atom)
		vals = [list(B_per_atom[at].values())[0] for at in ats]
		avg = np.average(vals)
		B_all = B_per_atom #{'0':avg}
		with open('struct%s_%s_%d_%s_rot.json' % (nstr,lstr,i_d,lint_str),'w') as writejson:
			json.dump(B_all,writejson, sort_keys=False, indent=2)
