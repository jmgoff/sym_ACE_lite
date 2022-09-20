from young import * 

# class for defining a coupling tree. It
#  contains tools to sort the tree within 
#  specified orbits of a permutation in S_N.
#  It is often referred to as a representation because,
#  the partitions of S_N are isomorphic to irreps in S_N.
#  the Tree_ID class adopts one of these partitions
class Tree_ID:
	def __init__(self,l,L,maplflag=False):
		self.l=l
		self.maplflag = maplflag
		if self.maplflag:
			self.map_l()
		self.L = L
		self.rank = len(l)
		self.set_sym_block()
		self.set_sigma_c_parts()
		self.orbit = None
		self.orbit_id = None

	def map_l(self):
		unique_ls = [u for u in sorted(list(set(self.l)))]
		lmap = {unique_ls[ind]:ind for ind in range(len(unique_ls))}
		revmapl = {ind:unique_ls[ind] for ind in range(len(unique_ls))}
		thisl = [lmap[i] for i in self.l]
		self.lmap = lmap
		self.revmapl = revmapl
		self.l = thisl

	def set_sym_block(self):
		# sets the highest symmetry partition for a given rank 
		#  it goes in powers of 2^x with some remainder. 
		sym_blocks=  {  4:((4,),),
			5:((4,1),),
			6:((4,2),),
			7:((4,2,1),),
			8:((8,),),
			}
		self.sym_block = sym_blocks[self.rank][0]

	def set_sigma_c_parts(self):
		ysgi = Young_Subgroup(self.rank)
		sigma_c_parts = ysgi.sigma_c_partitions(max_orbit=self.rank)
		sigma_c_parts.sort(key=lambda x: x.count(2),reverse=True)
		self.sigma_c = sigma_c_parts[0]
		sigma_c_parts.sort(key=lambda x: tuple([i%2==0 for i in x]),reverse=True)
		sigma_c_parts.sort(key=lambda x: max(x),reverse=True)
		self.sigma_c_parts = sigma_c_parts
		self.ysg = ysgi
		# generate leaf nodes structure for the binary tree
		nodes,remainder = tree(self.l)
		self.nodes = nodes
		self.remainder = remainder

	def return_leaf_l_ID(self):
		#if len(self.l) < 6:
		leaf_l = group_vec_by_node(self.l , self.nodes , self.remainder)
		leaf_id = [(tuple(sorted(lli)),Li) for lli,Li in zip(leaf_l,self.L)]
		leaf_only = [tuple(sorted(lli)) for lli in leaf_l]
		leaf_id = tuple(sorted(leaf_id))
		self.ltree = leaf_id
		self.leaf_only_l = leaf_only
		self.l_leaf_sym = [len(set(ln)) == 1 for ln in leaf_l]
		return leaf_id
	
	def return_orbit_l_ID(self,orbit=None):
		if orbit == None and self.orbit == None:
			orbit = self.sym_block
			self.orbit = orbit
		elif orbit != None:
			self.orbit = orbit
		else:
			orbit = self.orbit
		orb_l = group_vec_by_orbits(self.l, orbit)
		orb_id = []
		orb_only = []
		Lcount = 0
		for lorbit in orb_l:
			L_count_add = math.ceil(len(lorbit)/2)
			orbit_L = self.L[Lcount:Lcount + L_count_add] 
			Lcount += L_count_add
			orbit_nodes, orbit_remainder = tree(lorbit)
			leaf_l_perorb = group_vec_by_node(lorbit , orbit_nodes , orbit_remainder)
			orbit_leaf_id = tuple(sorted(leaf_l_perorb))
			orbit_leaf_id = sorted(orbit_leaf_id)
			orbit_leaf_add = tuple(sorted(lorbit))
			if self.maplflag:
				orbit_leaf_add = tuple([self.revmapl[ikl] for ikl in orbit_leaf_add])
			orb_id.append(tuple(orbit_leaf_add))
		orb_id = tuple(orb_id)
		self.orbit_id = orb_id
		self.l_orb_sym = [len(set(ln)) == 1 for ln in orb_l]
		return orb_id
	

	def return_orbit_nl_ID(self,nin,orbit =None):
		unique_ns = [u for u in sorted(list(set(nin)))]
		nmap = {unique_ns[ind]:ind for ind in range(len(unique_ns))}
		revmap = {ind:unique_ns[ind] for ind in range(len(unique_ns))}
		n = [nmap[i] for i in nin]
		if self.orbit == None:
			new_lid = self.return_orbit_l_ID(orbit=orbit)
		if orbit != None:
			self.orbit = orbit
			new_lid = self.return_orbit_l_ID(orbit=orbit)
		new_lid = self.orbit_id
		orb_n = group_vec_by_orbits(n, self.orbit)
		onl_id = []
		degen_id = []
		for ol,on in zip(new_lid,orb_n):
			itups = [(nii,lii) for nii,lii in zip(on,ol)]
			itups = sorted(itups)
			degen_id.append(tuple([(nii,lii) for nii,lii in zip(sorted(on),sorted(ol))]))
			orb_inner = tuple([ii[0] for ii in itups])
			id_inner = tuple([ii[1] for ii in itups])
			onl=tuple([orb_inner,id_inner])
			onl_id.append(onl)
		onl_id = sorted(onl_id)
		degen_id = tuple(sorted(degen_id))
		self.degen_id = degen_id
		onl_id.sort(key=lambda x: len(x),reverse = True)
		onl_id.sort(key=lambda x: len(x[0]),reverse = True)
		return tuple(onl_id)

