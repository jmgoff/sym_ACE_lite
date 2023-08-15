import os
os.environ['NUMPY_EXPERIMENTAL_ARRAY_FUNCTION'] = '0'
import scipy
from multiprocessing import pool,cpu_count
from gen_labels import *
from descriptor_calc_local.site_basis import *
from descriptor_calc_local.convert_configs import *
from wigner_couple import *
from scipy import special
from scipy.spatial.transform import Rotation as R
import matplotlib.pyplot as plt
from ase.neighborlist import primitive_neighbor_list
from ase import Atoms,Atom
from ase.io import read,write

def combine_ccs(ccs,linters):
	#global_ccs[rank][lstr][linter] - > localinput=ccs
	# ccs[linter][mstr]
	combined = { }
	for linter in ccs.keys():
		for mstr in ccs[linter].keys():
			if mstr not in combined.keys():
				combined[mstr] = 0.
			combined[mstr] += ccs[linter][mstr]
	return combined

class RI:
	def __init__(self,
			angular,
			lvec,
			**kwargs):
		self.ab = angular
		ls = [int(k) for k in lvec.split(',')]
		self.rank = len(ls)
		self.ls = ls 
		self.b = 0.
		self.ccs = None
		self.linter = None
		calc_dct = {}
		try:
			self.linter = kwargs['linter']
		except KeyError:
			pass
		try:
			self.set_ccs(kwargs['ccs'])
		except KeyError:
			pass
		if self.linter !=None:
			calc_dct['linter'] = self.linter
		if self.ccs !=None:
			calc_dct['ccs'] = self.ccs
		return None

	def set_ccs(self,ccs,linter=None):
		self.ccs = ccs

	def calc_invariant(self,**kwargs):
		if self.ccs !=None:
			ccs = self.ccs
		else:
			ccs = global_ccs
		if self.rank ==1:
			l = self.ls
			self.b = my_sum(self.ab.ylm(l,0)) #hard coded for M_R=0 for now
		elif self.rank >1:
			llst = ['%d']*self.rank
			lstr = ','.join(b for b in llst) % tuple(self.ls)
			if self.linter == None:
				mcombstrs = list(ccs[self.rank][lstr].keys())
			elif self.linter != None:
				print (self.ls,ccs[self.rank][lstr].keys())
				mcombstrs = list(ccs[self.rank][lstr][self.linter].keys())
			#sum over m combos
			A_precomp =  {l: {m:my_sum(self.ab.ylm(l,m)) for m in range(-l,l+1)} for l in range(max(self.ls)+1)} 
			for mcombstr in mcombstrs:
				mlst = [int(m) for m in mcombstr.split(',')]
				mlcheck = [mlst[i] in range(-self.ls[i],self.ls[i]+1) for i in range(len(mlst))]
				if all(mlcheck):
					A_lst = [A_precomp[self.ls[i]][mlst[i]] for i in range(len(self.ls))]
				elif not all(mlcheck):
					A_lst = [0.]*len(mlst)
				if self.linter == None:
					self.b += self_prd(A_lst) * ccs[self.rank][lstr][mcombstr]
				elif self.linter != None:
					self.b += self_prd(A_lst) * ccs[self.rank][lstr][self.linter][mcombstr]

def atoms_RI(args,**kwargs):
	i_d = args['id']
	atoms = args['atoms']
	elements = args['elements']
	print (i_d,atoms)
	rc = args['rc']
	nradbase=args['nradbase']
	nradmax_dict =args['nradmax']
	lmax_dict = args['lmax']
	lmbda=args['lmbda']
	ranks=args['rank']
	nus=args['nus']
	ccs = None
	linter = None
	eb_dct = {}
	if 'ccs' in kwargs.keys():
		ccs = kwargs['ccs']
		eb_dct['ccs'] = ccs
	if 'linter' in kwargs.keys():
		linter = kwargs['linter']
		eb_dct['linter']=linter

	single_elm = True
	try:
		single_elm = kwargs['single_element']
	except KeyError:
		pass
	tol = 0.0

	#Build the neighbor list for the atoms (including periodic images if applicable)
	nl = primitive_neighbor_list('ijdD',pbc=atoms.pbc,positions=atoms.positions ,cell=atoms.get_cell(),cutoff=rc+tol)
	elems = [e for e in atoms.symbols]
	unique_elems = list(set(elems))
	unique_elems = sorted(unique_elems)
	mump = {r:u for r,u in zip(range(len(unique_elems)),unique_elems)}
	mumpr = {u:r for r,u in zip(range(len(unique_elems)),unique_elems)}
	atinds = [atom.index for atom in atoms]
	at_neighs = { i: [] for i in atinds}
	at_dists = {i:[] for i in atinds}
	at_mus = {i:[] for i in atinds}
	for i,j in zip(nl[0],nl[1]):
		at_neighs[i].append(j)
		at_mus[i].append(atoms[j].symbol)
	for i,j in zip(nl[0],nl[-1]):
		at_dists[i].append(j)

	try:
		transform = kwargs['transformation']
		for i in [atom.index for atom in atoms]:
			positions = at_dists[i]
			positions = np.array(positions)
			pos_rot = np.matmul(transform.as_matrix(),positions.T)
			at_dists[i] =pos_rot.T
	except KeyError:
		at_dists = at_dists
	
	#nu_linters = nus_lints_from_nus(nu
	nu_linters = []
	for lstr in nus:
		l = [int(k) for k in lstr.split(',')]
		linters = tree_l_inters(l)
		print (l,linters)
		for linter in linters:
			lint_lst = ['%d'] * len(linter)
			lint_str = '-'.join(k for k in lint_lst)
			lint_str = lint_str % linter
			lL_str = lstr + '_' + lint_str
			nu_linters.append(lL_str)
	print (nu_linters)
	B_per_atom = {atid : {nu: 0. for nu in nu_linters}  for atid in range(len(atoms))}
	for nu in nu_linters:
		lstr = nu.split('_')[0]
		linterstr = nu.split('_')[-1]
		linter = tuple([int(b) for b in linterstr.split('-')])
		rank = len(l)
		for atind in at_dists.keys():
			r = np.array([np.linalg.norm(p) for p in at_dists[atind]])
			ab = angular_basis(at_dists[atind],lmax=lmax_dict[rank])
			eb_dct['linter'] = linter
			ebi= RI(ab,lstr, **eb_dct)
			ebi.calc_invariant( **eb_dct)
			B_per_atom[atind][nu] = np.real(ebi.b)
			del ebi
	return B_per_atom
