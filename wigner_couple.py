import numpy as np
from wigner_tree import *
from sympy.combinatorics import Permutation
from functools import partial

def get_wig_coupling(ldict,L_R=0,**kwargs):
	M_Rs = list(range(-L_R,L_R+1))
	#generic coupling for any L_R - support must be added to call 
	ranks = list(ldict.keys())
	coupling = {M_R : {rank:{} for rank in ranks} for M_R in M_Rs}

	for M_R in M_Rs:
		for rank in ranks:
			rnk = rank
			ls_per_rnk = generate_l_LR(range(ldict[rank]+1),rank,L_R,M_R)
			for lstr in ls_per_rnk:
				l = [int(k) for k in lstr.split(',')]
				if rank ==1:
					decomped = rank_1_tree(l,L_R,M_R)
					coupling[M_R][rnk][lstr] = decomped
				elif rank ==2:
					decomped = rank_2_tree(l,L_R,M_R)
					coupling[M_R][rnk][lstr] = decomped
				elif rank ==3:
					decomped = rank_3_tree(l,L_R,M_R)
					coupling[M_R][rnk][lstr] = decomped
				elif rank ==4:
					decomped = rank_4_tree(l,L_R,M_R)
					coupling[M_R][rnk][lstr] = decomped
				elif rank ==5:
					decomped = rank_5_tree(l,L_R,M_R)
					coupling[M_R][rnk][lstr] = decomped
				elif rank ==6:
					decomped = rank_6_tree(l,L_R,M_R)
					coupling[M_R][rnk][lstr] = decomped
				elif rank ==7:
					decomped = rank_7_tree(l,L_R,M_R)
					coupling[M_R][rnk][lstr] = decomped
				elif rank ==8:
					decomped = rank_8_tree(l,L_R,M_R)
					coupling[M_R][rnk][lstr] = decomped
				elif rank > 8:
					raise ValueError("Cannot generate couplings for rank %d. symmetric L_R couplings up to rank 6 have been implemented" % rank)
	return coupling

