import os
os.environ['NUMPY_EXPERIMENTAL_ARRAY_FUNCTION'] = '0'
import scipy
from rpi_lib import *
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

#energy invariant
class eB:
	def __init__(self,
			elements, # all possible elements (will be sorted alphabetically)
			radial,
			angular,
			rank,
			nu,
			**kwargs):
		self.rb = radial
		self.ab = angular
		self.rank = rank
		self.nu = nu
		mu0,mus,ns,ls = get_mu_n_l(nu)
		self.ns = ns #vector
		self.ls = ls #vector
		self.mump = {ind:val for ind,val in enumerate(sorted(elements))}
		self.mu=','.join(b for b in [self.mump[m] for m in mus]) # muj in the descriptor label vector
		self.B = 0.
		self.mu0 = mu0
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
		#self.calc_invariant(**calc_dct)
		return None

	def set_ccs(self,ccs,linter=None):
		self.ccs = ccs

	def calc_invariant(self,**kwargs):
		if self.ccs !=None:
			ccs = self.ccs
		else:
			ccs = global_ccs
		if self.rank ==1:
			n = self.ns[0]
			l = self.ls
			alst = [A_basis(self.rb,self.ab)]
			A_lst = [a.A_1(n, **{'muj':self.mu}) for a in alst]
			self.B = np.prod(A_lst,dtype=numpy.double)
		elif self.rank >1:
			try:
				mulst = [muj for muj in self.mu.split(',')]
			except AttributeError:
				mulst = self.mu
			A_bas = A_basis(self.rb,self.ab)
			alst = [A_bas]*self.rank
			llst = ['%d']*self.rank
			lstr = ','.join(b for b in llst) % tuple(self.ls)
			if self.linter == None:
				mcombstrs = list(ccs[self.rank][lstr].keys())
			elif self.linter != None:
				mcombstrs = list(ccs[self.rank][lstr][self.linter].keys())

			A_precomp = {mu: {n: {l: {m:0. for m in range(-l,l+1)} for l in range(max(self.ls)+1)} for n in range(1,max(self.ns)+1)} for mu in list(set(mulst))}
			for mu in list(set(mulst)):
				for n in range(1,max(self.ns)+1):
					for l in range(max(self.ls)+1):
						for m in range(-l,l+1):
							A_precomp[mu][n][l][m] = A_bas.A(n,l,m, **{'muj':mu})

			for mcombstr in mcombstrs:
				mlst = [int(m) for m in mcombstr.split(',')]
				mlcheck = [mlst[i] in range(-self.ls[i],self.ls[i]+1) for i in range(len(mlst))]
				if all(mlcheck):
					A_lst = [A_precomp[mulst[i]][self.ns[i]][self.ls[i]][mlst[i]] for i,a in enumerate(alst)]
				elif not all(mlcheck):
					A_lst = [0.]*len(mlst)
				if self.linter == None:
					self.B += np.prod(A_lst) * ccs[self.rank][lstr][mcombstr]
				elif self.linter != None:
					self.B += np.prod(A_lst) * ccs[self.rank][lstr][self.linter][mcombstr]

	def calc_equivariant(self,L_R,M_R,**kwargs):
		local_cc_flag = False
		try:
			self.linter = kwargs['linter']
		except KeyError:
			self.linter = None
		if L_R == 0:
			if self.ccs == None:
				ccs = global_ccs
			else:
				ccs = self.ccs
		elif L_R == 1:
			ccs = global_ccs_l1[M_R]
		if self.rank ==1:
			n = self.ns[0]
			l = L_R
			alst = [A_basis(self.rb,self.ab)]
			A_lst = [a.A(n,l,M_R, **{'muj':self.mu}) for a in alst]
			self.B = np.prod(A_lst)#np.prod(A_lst)
		elif self.rank >1:
			try:
				mulst = [muj for muj in self.mu.split(',')]
			except AttributeError:
				mulst = self.mu
			A_bas = A_basis(self.rb,self.ab)
			alst = [A_bas]*self.rank
			llst = ['%d']*self.rank
			lstr = ','.join(b for b in llst) % tuple(self.ls)
			if self.linter == None:
				try:
					mcombstrs = list(ccs[self.rank][lstr][()].keys())
					self.linter = ()
				except KeyError:
					mcombstrs = list(ccs[self.rank][lstr].keys())
			elif self.linter != None:
				try:
					mcombstrs = list(ccs[self.rank][lstr][self.linter].keys())
				except KeyError:
					try:
						mcombstrs = list(ccs[self.rank][lstr][self.linter[0]].keys())
						self.linter = self.linter[0]
					except KeyError:
						local_cc_flag = True
			#sum over m combos
			A_precomp = {mu: {n: {l: {m:0. for m in range(-l,l+1)} for l in range(max(self.ls)+1)} for n in range(1,max(self.ns)+1)} for mu in list(set(mulst))}
			for mu in list(set(mulst)):
				for n in range(1,max(self.ns)+1):
					for l in range(max(self.ls)+1):
						for m in range(-l,l+1):
							A_precomp[mu][n][l][m] = A_bas.A(n,l,m, **{'muj':mu})
			if not local_cc_flag:
				for mcombstr in mcombstrs:
					mlst = [int(m) for m in mcombstr.split(',')]
					A_lst = [A_precomp[mulst[i]][self.ns[i]][self.ls[i]][mlst[i]] for i in range(self.rank)]
					if self.linter == None:
						self.B += np.prod(A_lst) * ccs[self.rank][lstr][mcombstr]
					elif self.linter != None:
						self.B += np.prod(A_lst) * ccs[self.rank][lstr][self.linter][mcombstr]
			elif local_cc_flag:
				print ('WARNING: calculating coupling coefficients on-the-fly - this may be slow')
				assert L_R ==0, "on-the-fly coupling coefficients not implemented for L_R = 1"
				local_ccs_all,wts = local_coupling(self.ls)
				local_ccs = local_ccs_all[0] 
				mcombstrs = list(local_ccs[self.rank][lstr][self.linter].keys())
				for mcombstr in mcombstrs:
					mlst = [int(m) for m in mcombstr.split(',')]
					A_lst = [A_precomp[mulst[i]][self.ns[i]][self.ls[i]][mlst[i]] for i in range(self.rank)]
					if self.linter == None:
						self.B += np.prod(A_lst) * local_ccs[self.rank][lstr][mcombstr]
					elif self.linter != None:
						self.B += np.prod(A_lst) * local_ccs[self.rank][lstr][self.linter][mcombstr]

def atoms_eB_0rep(args,**kwargs):
	i_d = args['id']
	atoms = args['atoms']
	elements = args['elements']
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
	
	nu_linters = nus.copy()

	B_per_atom = {atid : {nu: 0. for nu in nu_linters}  for atid in range(len(atoms))}
	for nu in nus:
		mu0,muinds,n,l,Lv = get_mu_n_l(nu,return_L=True)
		treels = [ tuple(Lv) ]
		for linter in treels:
			nu_lint = nu
			if type(linter) == tuple:
				lint_lst = ['%d'] * len(linter)
				lint_str = '-'.join(k for k in lint_lst)
				lint_str = lint_str % linter
			elif type(linter) == int:
				lint_str = str(linter)
			if len(lint_str) >= 1:
				nu_lint = nu_lint + '_' + lint_str
			rank = get_mu_nu_rank(nu)
			muspec = [mump[m] for m in muinds]
			musp = ','.join(b for b in muspec)
			for mu in [musp]:
				for atind in at_dists.keys():
					if mu0 == mumpr[elems[atind]]:
						r = np.array([np.linalg.norm(p) for p in at_dists[atind]],dtype=np.float64)
						neigh_specs = [atoms[atind].symbol for atind in list(at_dists.keys())]
						rb = radial_basis(r_arr = r, rc=rc, nradbase=nradbase, nradmax=nradmax_dict[rank], lmax=lmax_dict[rank], lmbda=lmbda)
						rb.set_species(at_mus[atind])
						if rank ==1:
							rb = radial_basis(r_arr = r, rc=rc, nradbase=nradbase, nradmax=nradbase, lmax=lmax_dict[rank], lmbda=lmbda)
							rb.set_species(at_mus[atind])
						ab = angular_basis(at_dists[atind],lmax=lmax_dict[rank])
						eb_dct['linter'] = linter
						ebi= eB(elements,rb,ab,get_mu_nu_rank(nu),nu, **eb_dct)
						B_per_atom[atind][nu_lint] = np.real(ebi.B)
	#remove descriptor labels with nonmatching mu0s (these arent calculated but the labels are still in the dictionary)
	simplified_B_per_atom = {atind: {} for atind in atinds}
	for atind in atinds:
		tmpnus = list(B_per_atom[atind].keys())
		for nu in tmpnus:
			mu0 = int(nu.split('_')[0])
			if mu0 == mumpr[elems[atind]]:
				simplified_B_per_atom[atind][nu] = B_per_atom[atind][nu]
			else:
				pass
	return simplified_B_per_atom

def atoms_eB(args,**kwargs):
	manual_tree = False
	i_d = args['id']
	atoms = args['atoms']
	elements = args['elements']
	rc = args['rc']
	nradbase=args['nradbase']
	nradmax_dict =args['nradmax']
	lmax_dict = args['lmax']
	L_R = args['L_R']
	M_R = args['M_R']
	print (i_d,atoms,L_R,M_R)
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
	at_diststmp = {i:[] for i in atinds}
	at_dists = {i:[] for i in atinds}
	at_mus = {i:[] for i in atinds}
	for i,j in zip(nl[0],nl[1]):
		at_neighs[i].append(j)
		at_mus[i].append(atoms[j].symbol)
	for i,j in zip(nl[0],nl[-1]):
		at_diststmp[i].append(j)

	try:
		transform = kwargs['transformation']
		for i in [atom.index for atom in atoms]:
			positions = at_diststmp[i]
			positions = np.array(positions)
			pos_rot = np.matmul(transform.as_matrix(),positions.T)
			at_dists[i] =pos_rot.T
	except KeyError:
		at_dists = at_diststmp
	nu_linters = nus.copy()

	B_per_atom = {atid : {nu: 0. for nu in nu_linters}  for atid in range(len(atoms))}
	for nu in nus:
		if manual_tree:
			mu0,muinds,n,l,Lv = get_mu_n_l(nu,return_L=True)
			try:
				treels = [ tuple(Lv) ]
			except TypeError:
				treels = [ None ]
		elif not manual_tree:
			mu0,muinds,n,l = get_mu_n_l(nu,return_L=False)
			treels = tree_l_inters(l,L_R)
		for linter in treels:
			nu_lint = nu
			rank = get_mu_nu_rank(nu)
			muspec = [mump[m] for m in muinds]
			musp = ','.join(b for b in muspec)
			for mu in [musp]:
				for atind in at_dists.keys():
					if mu0 == mumpr[elems[atind]]:
						r = np.array([np.linalg.norm(p) for p in at_dists[atind]])
						neigh_specs = [atoms[atind].symbol for atind in list(at_dists.keys())]
						rb = radial_basis(r_arr = r, rc=rc, nradbase=nradbase, nradmax=nradmax_dict[rank], lmax=lmax_dict[rank], lmbda=lmbda)
						rb.set_species(at_mus[atind])
						if rank ==1:
							rb = radial_basis(r_arr = r, rc=rc, nradbase=nradbase, nradmax=nradbase, lmax=lmax_dict[rank], lmbda=lmbda)
							rb.set_species(at_mus[atind])
						ab = angular_basis(at_dists[atind],lmax=lmax_dict[rank])
						eb_dct['linter'] = linter
						ebi = eB(elements,rb,ab,get_mu_nu_rank(nu),nu, **eb_dct)
						ebi.calc_equivariant(L_R,M_R,**{'linter':linter})
						B_per_atom[atind][nu_lint] = np.float64(ebi.B)
	#remove descriptor labels with nonmatching mu0s (these arent calculated but the labels are still in the dictionary)
	simplified_B_per_atom = {atind: {} for atind in atinds}
	for atind in atinds:
		tmpnus = list(B_per_atom[atind].keys())
		for nu in tmpnus:
			mu0 = int(nu.split('_')[0])
			if mu0 == mumpr[elems[atind]]:
				simplified_B_per_atom[atind][nu] = B_per_atom[atind][nu]
			else:
				pass
	return simplified_B_per_atom
