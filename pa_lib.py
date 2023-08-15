from inter_set import *
from symmetric_grp_manip import *

def tree_labels(nin,lin,L_R=0,M_R=0):
    rank = len(lin)
    ysgi = Young_Subgroup(rank)
    if type(lin) != list:
        lin = list(lin)
    if type(nin) != list:
        nin = list(nin)
    # get possible unique l permutations based on degeneracy and coupling tree structure
    ysgi.subgroup_fill(lin, partitions=[local_sigma_c_partitions[rank][-1]], max_orbit = len(local_sigma_c_partitions[rank][-1]), sigma_c_symmetric=False, semistandard=False)
    lperms = ysgi.fills.copy()
    lperms = leaf_filter(lperms)
    if rank not in [4,8,16,32]:
        lperms_tmp = []
        used_hrep = []
        for lperm in lperms:
            hrep = get_highest_coupling_representation(tuple(lperm),tuple(lperms[0]))
            if hrep not in used_hrep:
                used_hrep.append(hrep)
                lperms_tmp.append(lperm)
            else:
                pass
                #print('omitting',lperm)
        lperms = lperms_tmp
        #print ('lperms after subtree filter',lperms)
    original_joint_span = {lp:[] for lp in lperms}
    orb_nls = []

    ls = lperms.copy()
    nps_per_l = {}

    # get n permutations per l permutation
    # this could equivalently be done with a search over S_N
    for lp in ls:
        rank = len(lp)
        original_span_SO3 = tree_l_inters(lp) #RI basis size
        degen_orbit, orbit_inds = get_degen_orb(lp) #PI basis size
        ysgi.subgroup_fill(nin,[degen_orbit],sigma_c_symmetric=False,semistandard=False)
        degen_fills = ysgi.fills.copy()
        maxdegen = max([len(ois) for ois in orbit_inds])
        sequential_degen_orbit, orbit_inds_s = enforce_sorted_orbit(orbit_inds)
        #if rank > 4 and maxdegen > math.ceil(rank/2):
        #    ysgi.subgroup_fill(nin,[degen_orbit],sigma_c_symmetric=False,semistandard=False)
        #else:
        ysgi.subgroup_fill(nin,[sequential_degen_orbit],sigma_c_symmetric=False,semistandard=False)
        nps_per_l[lp] = ysgi.fills.copy()
        original_joint_span[lp] = [(prd[0],lp,prd[1]) for prd in itertools.product(degen_fills,original_span_SO3)]

    nlabs = 0
    labels_per_lperm = {}
    #build all labels (unsorted trees)
    for l in ls:
        #print ('in l loop',l)
        l = list(l)
        subblock = []
        rank = len(l)
        inters = tree_l_inters(list(l),L_R=L_R,M_R=M_R)
        nperms = nps_per_l[tuple(l)]
        muperms = [tuple([0]*rank)]
        for inter in inters:
            if rank <= 5:
                if np.sum([inter[0]] + l[:2]) %2 ==0 and np.sum([inter[1]] + l[2:4]) %2 ==0:
                    for muperm in muperms:
                        for nperm in nperms:
                            if rank == 5:
                                orb_nls.append("0_%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d_%d-%d-%d"% (muperm + nperm+tuple(l) + inter))
                                subblock.append("0_%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d_%d-%d-%d"% (muperm + nperm+tuple(l) + inter))
                            elif rank ==4:
                                orb_nls.append("0_%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d_%d-%d"% (muperm + nperm+tuple(l) + inter))
                                subblock.append("0_%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d_%d-%d"% (muperm + nperm+tuple(l) + inter))
                            nlabs +=1
        labels_per_lperm[tuple(l)] = subblock


    block_sizes = {key:len(val) for key,val in labels_per_lperm.items()}
    all_labs = []
    labels_per_block = {block:[] for block in sorted(list(block_sizes.keys()))}
    counts_per_block = {block:0 for block in sorted(list(block_sizes.keys()))}
    # collect sorted trees only
    for block,labs in labels_per_lperm.items():
        used_ns = []
        used_ids = []
        for nu in labs:
            mu0,mutst,ntst,ltst, L = get_mu_n_l(nu,return_L=True)
            ltree = [(li,ni) for ni,li in zip(ntst,ltst)] #sort first on n
            tree_i =  build_tree(ltree,L,L_R)
            ttup = tree_i.full_tup()
            tid = tree_i.tree_id
            conds = tid not in used_ids # sorting is ensured in construction of trees
            if conds:
                if tuple(ntst) not in used_ns:
                    used_ns.append(tuple(ntst))
                used_ids.append(tid)
                labels_per_block[block].append(nu)
                counts_per_block[block] += 1
                all_labs.append(nu)
            else:
                pass

    #collect labels per l permutation block
    max_labs = []
    max_count = max(list(counts_per_block.values()))
    for block,tree_labs in labels_per_block.items():
        if len(tree_labs) == max_count:
            max_labs.append(tree_labs.copy())
    max_labs = max_labs[0]

    return max_labs,all_labs,labels_per_block,original_joint_span

def combine_blocks(blocks,lin,original_spans,L_R=0):
    # tool to recombine trees from multiple permutations of l
    rank = len(lin)
    ysgi = Young_Subgroup(rank)
    lps = list(blocks.keys())
    blockpairs = [(block1,block2) for block1,block2 in itertools.product(lps,lps) if block1 != block2]
    if len(blockpairs) == 0:
        blockpairs = [(block1,block2) for block1,block2 in itertools.product(lps,lps)]
    block_map = {blockpair:None for blockpair in blockpairs}
    all_map = {blockpair:None for blockpair in blockpairs}
    L_map = {blockpair:None for blockpair in blockpairs}
    raw_perms = [p for p in itertools.permutations(list(range(rank)))]
    Ps = [Permutation(pi) for pi in raw_perms]
    for blockpair in list(block_map.keys()):
        l1i,l2i = blockpair
        is_sigma0 = l1i == lps[0]
        Pl1is = [P for P in Ps if tuple(P(list(l1i))) == l2i]
        Pl1_maxorbit_sizes = [max([len(k) for k in P.full_cyclic_form]) for P in Pl1is]
        maxorbit_all = max(Pl1_maxorbit_sizes)
        maxorbit_ind = Pl1_maxorbit_sizes.index(maxorbit_all)
        if not is_sigma0:
            block_map[blockpair] = Pl1is[maxorbit_ind]
            all_map[blockpair] = Pl1is
        else:
            block_map[blockpair] = Permutation(tuple( [ tuple([ii]) for ii in list(range(rank))]))
            all_map[blockpair] = [Permutation(tuple( [ tuple([ii]) for ii in list(range(rank))]))]
            

    for blockpair in list(block_map.keys()):
        l1i,l2i = blockpair
        inters1 = tree_l_inters(l1i,L_R)
        is_sigma0 = tuple(l1i) == lps[0]
        l1i = list(l1i)
        l2i = list(l2i)
        # intermediates hard coded for ramk 4 and 5 right now
        inters1 = [inter for inter in inters1 if np.sum([inter[0]] + l1i[:2]) %2 ==0 and np.sum([inter[1]] + l1i[2:4]) %2 ==0 ]
        inters2 = tree_l_inters(l2i,L_R)
        inters2 = [inter for inter in inters2 if np.sum([inter[0]] + l2i[:2]) %2 ==0 and np.sum([inter[1]] + l2i[2:4]) %2 ==0 ]
        if not is_sigma0:
            L_map[blockpair] = {L1i:L2i for L1i,L2i in zip(inters1,inters2)  }
        else:
            L_map[blockpair] = {L1i:L1i for L1i,L1i in zip(inters1,inters1)  }
    used_ids = []
    used_nl = []
    combined_labs = []
    super_inters_per_nl = {}
    for lp,nus in blocks.items():
        rank = len(lp)
        degen_orbit, orbit_inds = get_degen_orb(lp)
        maxdegen = max([len(ois) for ois in orbit_inds])
        sequential_degen_orbit, orbit_inds_s = enforce_sorted_orbit(orbit_inds)
        block_pairs = [blockpair for blockpair in list(block_map.keys()) if blockpair[0] == tuple(lp)]
        blockpair = block_pairs[0]
        #perm_map = block_map[block_pairs[0]]
        if rank == 4:
            perms_2_check = [block_map[blockpair]]
        else:
            perms_2_check = [block_map[blockpair]]
            #perms_2_check = all_map[blockpair]
        for nu in nus:
            mu0ii,muii,nii,lii,Lii = get_mu_n_l(nu,return_L=True)
            is_sigma0 = tuple(lii) == lps[0]
            degen_orbit, orbit_inds = get_degen_orb(lp)
            sequential_degen_orbit, orbit_inds_s = enforce_sorted_orbit(orbit_inds)
            nlii = [(niii,liii) for niii,liii in zip (nii,lii)]
            atrees = []
            for perm_map in perms_2_check:
                remapped = perm_map(nlii)
                newnii = [nliii[0] for nliii in remapped] 
                newlii = [nliii[1] for nliii in remapped]
                new_Lii = L_map[blockpair][Lii]
                new_ltree = [(liii,niii) for niii,liii in zip(newnii,newlii)] 
                tree_i =  build_tree(new_ltree,Lii,L_R)
                #tree_i =  build_tree(new_ltree,new_Lii,L_R)
                ttup = tree_i.full_tup()
                tid = tree_i.tree_id
                atrees.append(tid)
            cond1 = not any([tid in used_ids for tid in atrees])
            #cond1 = tid not in used_ids
            if is_sigma0:
                #testing
                #cond2 = (tuple(newnii),tuple(newlii)) not in used_nl
                cond2 = True 
            else:
                #cond2 = (tuple(newnii),tuple(newlii)) not in used_nl        
                cond2 = True 

            if cond1 and cond2:
                combined_labs.append(nu)
                used_ids.append(tid)
                used_nl.append((tuple(newnii),tuple(newlii)))
                try:
                    super_inters_per_nl[(tuple(newnii),tuple(newlii))].append(new_Lii)
                except KeyError:
                    super_inters_per_nl[(tuple(newnii),tuple(newlii))] = [new_Lii]
            else:
                pass
    return combined_labs


seq_degen_map = {
((2,2,1),(5,)):(4,1),

((2,2,1),(4,1)):(4,1),
((2,1,1,1),(4,1)):(3,1,1),

((2,2,1),(2,3)):(2,2,1),
((2,1,1,1),(2,3)):(2,1,1,1),
((1,1,1,1,1),(2,3)):(2,1,1,1),
}

#apply ladder relationships
def apply_ladder_relationships(lin, nin, combined_labs, parity_span, parity_span_labs, full_span, L_R=0):
    N = len(lin)
    uniques = list(set(lin))
    tmp = list(lin).copy()
    tmp.sort(key=Counter(lin).get,reverse=True)
    uniques.sort()
    uniques.sort(key=Counter(tmp).get,reverse=True)
    count_uniques =[lin.count(u) for u in uniques]
    mp = {uniques[i]:i for i in range(len(uniques))}
    mprev = {i:uniques[i] for i in range(len(uniques))}
    mappedl = [mp[t] for t in tmp]
    ysgi = Young_Subgroup(N)

    unique_ns =  list(set(nin))
    tmpn = list(nin).copy()
    tmpn.sort(key=Counter(nin).get,reverse=True)
    unique_ns.sort()
    unique_ns.sort(key=Counter(nin).get,reverse=True)
    count_unique_ns =[nin.count(u) for u in unique_ns]
    mp_n = {unique_ns[i]:i for i in range(len(unique_ns))}
    mprev_n = {i:unique_ns[i] for i in range(len(unique_ns))}
    mappedn = [mp_n[t] for t in tmpn]
    mappedn = tuple(mappedn)
    mappedl = tuple(mappedl)

    max_labs = parity_span_labs.copy()
    #mapldegenrep, maplorbit_inds = get_degen_orb(mappedl)
    #ndegen_rep = list(ndegen_rep)
    #ndegen_rep.sort(key=lambda x: x, reverse =True)
    #ndegen_rep = tuple(ndegen_rep)
    # get partition of S_N that the vector of n are commensurate with
    #  based on degeneracy
    ndegen_rep, n_orbit_inds = get_degen_orb(mappedn)
    origndegen_rep, orign_orbit_inds = get_degen_orb(nin)
    ndegen_rep = list(ndegen_rep)
    ndegen_rep.sort(key=lambda x: x, reverse =True)
    ndegen_rep = tuple(ndegen_rep)
    degen_fam = (mappedl,ndegen_rep)

    all_inters = tree_l_inters(lin)
    even_inters = simple_parity_filt(lin, all_inters, L_R)

    if 0 in lin:
        funcs = combined_labs[:len(full_span)]

    else:
        if degen_fam ==   ((0,0,0,0),(4,)):
            funcs = parity_span_labs[::3] #according to full degen ladder relationship
        elif degen_fam == ((0,0,0,0),(3,1)):
            funcs = parity_span_labs[::3]
        elif degen_fam == ((0,0,0,0),(2,2)):
            funcs = parity_span_labs[:len(parity_span)]
        elif degen_fam == ((0,0,0,0),(2,1,1)):
            funcs = parity_span_labs[:len(parity_span)]
        elif degen_fam == ((0,0,0,0),(1,1,1,1)):
            funcs = combined_labs[:len(full_span)]

        elif degen_fam == ((0,0,0,1),(4,)):
            funcs = []
            recurmax = len(max_labs)/2
            count = 0
            for lab in max_labs:
                mu0ii, muii, nii, lii, Lii = get_mu_n_l(lab,return_L=True)
                lidegen_rep, l_orbit_inds = get_degen_orb(lii)
                ysgi.subgroup_fill(list(nin),[lidegen_rep],sigma_c_symmetric=False,semistandard=False)
                degen_nfills = ysgi.fills.copy()
                if count < recurmax and tuple(nii) in degen_nfills:
                    funcs.append(lab)
                    count += 1
        elif degen_fam == ((0,0,0,1),(3,1)):
            funcs = []
            recurmax = len(max_labs)/2
            count = 0
            for lab in max_labs:
                mu0ii, muii, nii, lii, Lii = get_mu_n_l(lab,return_L=True)
                lidegen_rep, l_orbit_inds = get_degen_orb(lii)
                ysgi.subgroup_fill(list(nin),[lidegen_rep],sigma_c_symmetric=False,semistandard=False)
                degen_nfills = ysgi.fills.copy()
                if count < recurmax and tuple(nii) in degen_nfills:
                    funcs.append(lab)
                    count += 1
        elif degen_fam == ((0,0,0,1),(2,2)):
            funcs = parity_span_labs[:len(parity_span)]
        elif degen_fam == ((0,0,0,1),(2,1,1)):
            funcs = []
            recurmax = len(max_labs)/2
            count = 0
            for lab in max_labs:
                mu0ii, muii, nii, lii, Lii = get_mu_n_l(lab,return_L=True)
                lidegen_rep, l_orbit_inds = get_degen_orb(lii)
                l_sequential_degen_orbit, l_orbit_inds_s = enforce_sorted_orbit(l_orbit_inds)
                # switch to lower symmetry SN representation
                ysgi.subgroup_fill(list(nin),[l_sequential_degen_orbit],sigma_c_symmetric=False,semistandard=False)
                degen_nfills = ysgi.fills.copy()
                if count < recurmax and tuple(nii) in degen_nfills:
                    funcs.append(lab)
                    count += 1
        elif degen_fam == ((0,0,0,1),(1,1,1,1)):
            funcs = combined_labs[:len(full_span)]

        elif degen_fam == ((0,0,1,1),(4,)):
            funcs = parity_span_labs
        elif degen_fam == ((0,0,1,1),(3,1)):
            funcs = parity_span_labs
        elif degen_fam == ((0,0,1,1),(2,2)):
            funcs = combined_labs[:len(parity_span) + len(even_inters[1:])]
        elif degen_fam == ((0,0,1,1),(2,1,1)):
            funcs = combined_labs[:len(parity_span) + (2*len(even_inters[1:]))]
        elif degen_fam == ((0,0,1,1),(1,1,1,1)):
            funcs = combined_labs[:len(full_span)]

        elif degen_fam == ((0,0,1,2),(4,)):
            funcs = parity_span_labs
        elif degen_fam == ((0,0,1,2),(3,1)):
            funcs = combined_labs[:len(parity_span) + len(even_inters[1:])]
        elif degen_fam == ((0,0,1,2),(2,2)):
            funcs = combined_labs[:len(parity_span) + len(all_inters[1:])]
        elif degen_fam == ((0,0,1,2),(2,1,1)):
            funcs = combined_labs[:len(parity_span) + ((len(all_inters) -1)*2) + len(even_inters[1:])]
        elif degen_fam == ((0,0,1,2),(1,1,1,1)):
            funcs = combined_labs[:len(full_span)]

        elif degen_fam[0] == (0,1,2,3):
            funcs = combined_labs[:len(full_span)]



        elif degen_fam == ((0,0,0,0,0),(5,)):
            #print (len(full_span),len(parity_span),len(max_labs),len(even_inters))
            combined_labs.reverse()
            funcs = sorted(combined_labs[::4]) # from rank 5 ladder relationship

    return funcs

