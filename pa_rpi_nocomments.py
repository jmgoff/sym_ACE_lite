from lib.coupling_tree import * 
import sys


# code without comments to show simplicity.
debug = False

global_lsyms = {}
global_lsemi = {}
global_parts = {}
global_orbits = {}


def permutation_adapted_lL(l,semistandardflag=True):
	N = len(l)
	uniques = list(set(l))
	tmp = l.copy()
	tmp.sort(key=Counter(l).get,reverse=True)
	uniques.sort(key=Counter(tmp).get,reverse=True)
	count_uniques =[l.count(u) for u in uniques]
	mp = {uniques[i]:i for i in range(len(uniques))}
	mprev = {i:uniques[i] for i in range(len(uniques))}
	mapped = [mp[t] for t in tmp]
	
	if len(uniques) == math.floor(len(l)/2):# and 0 in uniques:
		l.sort(key=Counter(l).get,reverse=True)
		mapped.sort(key=Counter(mapped).get,reverse=True)
	else:
		l = sort_pair(l)
		mapped = sort_pair(mapped)
	assert len(l) <= 6, "symmetry reduction is only needed for rank 4 + descriptors. Use a rank between 4 and 6, automorphism groups for rank >= 7 to be added soon."
	N = len(l)
	if debug:
		print ('lrep',l)
	
	try:
		sigma_c_parts = global_parts[N]
	except KeyError:
		
		ysgi = Young_Subgroup(N)
		sigma_c_parts = ysgi.sigma_c_partitions(max_orbit=N)
		sigma_c_parts.sort(key=lambda x: x.count(2),reverse=True)
		sigma_c_parts.sort(key=lambda x: tuple([i%2==0 for i in x]),reverse=True)
		sigma_c_parts.sort(key=lambda x: max(x),reverse=True)
		global_parts[N] = sigma_c_parts
		if debug:
			print (sigma_c_parts)
	try:
		lperms = global_lsyms[tuple(mapped)]
		part_per_fill = global_orbits[tuple(mapped)]
		lperms_semistandard = global_lsemi[tuple(mapped)]
	except KeyError:
		ysgi = Young_Subgroup(N)
		ysgi.subgroup_fill(mapped,sigma_c_parts,sigma_c_symmetric=True,semistandard=False)
		lperms = ysgi.fills.copy()
		lperms = ysgi.reduce_list(mapped,lperms)
		part_per_fill = ysgi.partitions_per_fill
		ysgi.subgroup_fill(mapped,sigma_c_parts,sigma_c_symmetric=True,semistandard=True)
		lperms_semistandard = ysgi.fills.copy()
		standard_part_per_fill = ysgi.partitions_per_fill
		global_lsyms[tuple(mapped)] = lperms
		global_orbits[tuple(mapped)] = part_per_fill
		global_lsemi[tuple(mapped)] = lperms_semistandard
		if debug:
			print ('lperms_ysg_not',lperms)
			print ('lperms_ysg_semistandard',lperms_semistandard)
			print ('highest symmetry orbits per filling', part_per_fill)
	myvarsigma_l = [tuple([mprev[k] for k in lp]) for lp in lperms]
	myvarsigma_l = sorted(myvarsigma_l)
	varsigma_l = []

	L_per_varsigma = {}

	for vsl in myvarsigma_l:
		varsigma_l.append(vsl)
		L_lp = tree_l_inters(vsl)
		nodes,remainder = tree(vsl)
		lnodes_lp = group_vec_by_node(vsl,nodes,remainder)
		L_lp_filtered = parity_filter(lnodes_lp,L_lp)
		L_per_varsigma[ vsl]=L_lp_filtered

	reps_per_varsigma = {}
	for fill,irreps in part_per_fill.items():
		lmp = tuple([mprev[k] for k in fill])
		reps_per_varsigma[lmp] = irreps
	
	quick_build = []
	orbids = []
	leafids = []
	used_l = []
	all_autos = []
	ysg = Young_Subgroup(N)
	for lp, Lp_lst in L_per_varsigma.items():
		ti = Tree_ID(lp,Lp_lst[0])
		orbi = ti.return_orbit_l_ID()
		leafi = ti.return_leaf_l_ID()
		ysg.set_inds(lp)
		automorphisms = get_auto_part(tuple(lp),tuple([len(lp)]),add_degen_autos=False,part_only=False)
		conj1, applied_conj,conj_list = ysg.apply_automorphism_conjugation(my_automorphisms=automorphisms)
		if leafi not in leafids:
			quick_build.append( (lp,Lp_lst[0]))
			orbids.append(orbi)
			leafids.append(leafi)
			all_autos.extend(conj_list)
		if orbi not in orbids:
			quick_build.append( (lp,Lp_lst[0]))
			orbids.append(orbi)
			leafids.append(leafi)
			all_autos.extend(conj_list)

	varsigma_lL = quick_build.copy()
	highest_sym_reps = {}
	for lL in varsigma_lL:
		varsigmali,intersi = lL
		reps = reps_per_varsigma[varsigmali]
		highest_sym_reps[lL] = reps[0]
			
	return varsigma_lL, highest_sym_reps


def permutation_adapted_nlL(n,l,semistandardflag=True):
	return_desclabels = True
	lLs, SN_irreps = permutation_adapted_lL(l,semistandardflag)
	varsigma_nls = []
	used_nl_reps = []
	ns = unique_perms(n)
	
	for lL,irrep in SN_irreps.items():
		li,Li = lL
		ti = Tree_ID(li,Li)
		loid = ti.return_orbit_l_ID(orbit=irrep)
		for nperm in ns:
			nloid = ti.return_orbit_nl_ID(nperm,orbit=irrep)
			if debug:
				print (nloid)
			if nloid not in used_nl_reps:
				used_nl_reps.append(nloid)
				varsigma_nls.append((nperm,li,Li))

	

	descriptor_labels = []
	nelements=1
	N = len(l)
	for mu0 in range(nelements):
		for munlL in varsigma_nls:
			munlL = tuple([ tuple([0]*N)  ]) + munlL
			st='%d_' % mu0
			tmp=  ','.join([b for b in ['%d']*N*(3)])
			tmp = tmp % tuple(flatten(munlL[:-1]))
			st +=tmp
			st+=  '_'
			st+= '-'.join(str(b) for b in munlL[3])
			descriptor_labels.append(st)
	if return_desclabels:
		return descriptor_labels
	else:
		return varsigma_nls

def descriptor_labels_YSG(rank,nmax,lmax,mumax=0,lmin=1):
	if rank >= 4:
		lrng = list(range(lmin,lmax+1))
		nrng = list(range(1,nmax+1))
		lstrs = generate_l_LR(lrng , rank , L_R = 0)

		ns = [i for i in itertools.combinations_with_replacement(nrng,rank)]
		used_ls = []
		labels = []
		
		for lstr in sorted(list(set(lstrs))):
			l = [int(k) for k in lstr.split(',')]
			N = len(l)
			for n in ns:
				lrep = tuple([ tuple(sorted(l)), tuple(sorted(n))])
				if lrep not in used_ls:
					nls = permutation_adapted_nlL(n,l)
					used_ls.append(lrep)
					labels.extend(nls)
	elif rank < 4:
		labels = generate_nl(rank,nmax,lmax,mumax=mumax,lmin=lmin,L_R=0,M_R=0,all_perms=False)	
	return labels
