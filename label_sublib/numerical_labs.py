from gen_labels import *
from inter_set import *
from label_sublib.young import *
from rank5_PA_method import *

global_lex_labels_4 = {}

def lex_nl_labels_4(n,l,L_R=0):
        N = len(l)
        uniques = list(set(l))
        tmp = l.copy()
        tmp.sort(key=Counter(l).get,reverse=True)
        uniques.sort()
        uniques.sort(key=Counter(tmp).get,reverse=True)
        count_uniques =[l.count(u) for u in uniques]
        mp = {uniques[i]:i for i in range(len(uniques))}
        mprev = {i:uniques[i] for i in range(len(uniques))}
        mappedl = [mp[t] for t in tmp]

        unique_ns =  list(set(n))
        tmpn = n.copy()
        tmpn.sort(key=Counter(n).get,reverse=True)
        unique_ns.sort()
        unique_ns.sort(key=Counter(n).get,reverse=True)
        count_unique_ns =[n.count(u) for u in unique_ns]
        mp_n = {unique_ns[i]:i for i in range(len(unique_ns))}
        mprev_n = {i:unique_ns[i] for i in range(len(unique_ns))}
        mappedn = [mp_n[t] for t in tmpn]
        mappedn = tuple(mappedn)
        mappedl = tuple(mappedl)
        inters = tree_l_inters(tmp)
        #inters = simple_parity_filt(tmp, inters, L_R, even = True)
        lp_map = { (0,0,0,0):[(0,0,0,0)], #from 2,2 YSG
                (0,0,0,1):[(0,0,0,1)],
                (0,0,1,1):[(0,0,1,1),(0,1,0,1)],
                #(1,1,0,0):[(0,0,1,1)],
                (0,0,1,2):[(0,0,1,2),(0,1,0,2)],
                #(0,0,2,1):[(0,0,1,2)],
                (0,1,2,3):[(0,1,2,3),(0,2,1,3)],
                #(0,1,3,2):[(0,1,2,3)],
                #(0,3,1,2):[(0,1,2,3)],
                #(0,2,1,3):[(0,1,2,3)],
                #(0,3,2,1):[(0,1,2,3)],
                #(0,2,3,1):[(0,1,2,3)],
                #(1,3,2,0):[(0,1,2,3)],
                #(1,2,3,0):[(0,1,2,3)],
                #(2,3,0,1):[(0,1,2,3)],
                #(3,2,0,1):[(0,1,2,3)],
                #(3,2,1,0):[(0,1,2,3)],
                #(2,3,1,0):[(0,1,2,3)],
                #(1,2,0,3):[(0,1,2,3)],
                #(1,3,0,2):[(0,1,2,3)],
                #(1,3,2,0):[(0,1,2,3)],
                #(3,1,0,2):[(0,1,2,3)],
                #(3,1,2,0):[(0,1,2,3)],

                #(0,0,1,1):[(0,0,1,1),(0,1,0,1)],
                #(0,0,1,2):[(0,0,1,2),(0,1,0,2)],
                #(0,1,2,3):[(0,1,2,3),(0,2,1,3),(0,3,1,2)]
                }
        def procedure():
            #lps = lp_map[mappedl]
            lps = [mappedl]
            npstrs = ind_vec(mappedn,size=len(l))
            nps = [tuple([int(b) for b in npp.split(',')]) for npp in npstrs]
            nlps_all = [p for p in itertools.product(nps,lps)]
            nls_lex = []
            for nlp in nlps_all:
                npi,lpi = nlp
                tupped = tuple([(ni,li) for ni,li in zip(npi,lpi)])
                if tupped == tuple(sorted(tupped)):
                    nls_lex.append(nlp)

            nlLunmap = [p for p in itertools.product(nls_lex,inters)] 

            nlLs_partial = []
            for nlLu in nlLunmap:
                nli,Li  =nlLu
                lmp = [mprev[li] for li in nli[1]]
                nlLi_partial = (nli[0],tuple(lmp),Li)
                nlLs_partial.append(nlLi_partial)
            return nlLs_partial

        try:
            nlLs_partial = global_lex_labels_4[(mappedn,tuple(tmp))]
        except KeyError:
            nlLs_partial = procedure()
            global_lex_labels_4[(mappedn,tuple(tmp))] = nlLs_partial
        nlLs = []
        for nlLup in nlLs_partial:
            ni,li,Li  =nlLup
            nmp = [mprev_n[nii] for nii in ni]
            lmp = li
            nlLi = (tuple(nmp),lmp,Li)
            nlLs.append(nlLi)

        #print (nlLs)
        nus = []
        for nlL in nlLs:
            muvec = [0] * len(l)
            mustrlst = ['%d'] * len(l)
            nstrlst = ['%d'] * len(l)
            lstrlst = ['%d'] * len(l)
            mustr = ','.join(b for b in mustrlst) % tuple(muvec)
            np,lp,Lp = nlL
            Lstrlst = ['%d'] * len(Lp)
            Lstr = '-'.join(b for b in Lstrlst) % Lp
            mustr = ','.join(b for b in mustrlst) % tuple(muvec)
            nstr = ','.join(b for b in mustrlst) % np
            lstr = ','.join(b for b in mustrlst) % lp
            nu ='0_'+ mustr + ',' + nstr + ',' + lstr + '_' + Lstr
            nus.append(nu)
        return nus

def lexico_labels_nl(rank,nmax,lmax,L_R=0,M_R=0,lmin=0,return_SN=False):
        #if lmin > 0:
        #    print ('!!!! WARNING !!!! do not use this unless you will build products of ACE descriptors in your basis')

        assert L_R %2 ==0 ,'odd L_R not implemented yet.'
        nrng = range(1,nmax+1)
        nstrs = ind_vec(nrng,size=rank)
        lrng = range(lmin,lmax+1)
        nvecs = [tuple([int(b) for b in ni.split(',')]) for ni in nstrs]
        lstrs = generate_l_LR(lrng,rank,L_R,M_R)
        lvecs = [tuple([int(b) for b in li.split(',')]) for li in lstrs]
        used_nls = []
        nus = []
        all_SN = []
        #prds = []
        #for prd in itertools.product(nvecs,lvecs):
        #    tupped = tuple([(li,ni) for li,ni in zip(prd[0],prd[1])])
        #    if tupped == tuple(sorted(tupped)):
        #        prds.append(prd)
        prds = [p for p in itertools.product(nvecs,lvecs) if p[0] == tuple(sorted(p[0])) and p[1] == tuple(sorted(p[1]))]
        for prd in prds:
            n,l = prd
            if l != tuple(sorted(l)):
                raise ValueError(' problem with l %d,%d,%d,%d' % l)
            l = list(l)
            n = list(n)
            max_labs,all_labs,labels_per_block,original_spans = tree_labels(n,l)
            needed_lperms = list(labels_per_block.keys())
            num_needed_l = len(needed_lperms)
            vecstrlst = ['%d']*rank
            interstrlst = ['%d']*(rank-2)
            vecstr = ','.join(vecstrlst)
            interstr = '-'.join(interstrlst)
            for lperm in needed_lperms:
                lpl = list(lperm)
                these_ints = tree_l_inters(lpl,L_R)
                current_labels = lex_nl_labels_4(n,lpl,L_R)
                nus.extend(current_labels)
                sn_prds = [p for p in itertools.product(nvecs,[lperm],these_ints)]
                allsnprds = []
                for snprd in sn_prds:
                    nii,lii,Lii = snprd
                    mustrii = vecstr % tuple([0]*rank)
                    nstrii = vecstr % tuple(nii)
                    lstrii = vecstr % tuple(lii)
                    Lstrii = interstr % tuple(Lii)
                    nu = '0_' + mustrii + ',' + nstrii + ',' + lstrii + '_' + Lstrii
                    allsnprds.append(nu)
                    
                all_SN.extend(allsnprds)
            #print (prd, len(current_labels))
            #used_nls.append(prd)

        if return_SN:
            return nus,all_SN
        else:
            return nus


#tst = lex_nl_labels_4([1,1,1,2],[1,1,2,2],0)
#print (tst)
#lex_labs = lexico_labels_nl(rank=4,nmax=4,lmax=4,L_R=0,M_R=0,lmin=1)
#print (len(lex_labs))
