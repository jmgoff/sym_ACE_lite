from scipy import special
from sym_ACE_settings import *
import pickle
import numpy as np
import os


#lib_path = os.getcwd()

def Clebsch_gordan(j1,m1,j2,m2,j3,m3):
	# Clebsch-gordan coefficient calculator based on eqs. 4-5 of:
	# https://hal.inria.fr/hal-01851097/document

	#VERIFIED: test non-zero indices in Wolfram using format ClebschGordan[{j1,m1},{j2,m2},{j3,m3}]
	#rules:
	rule1 = np.abs(j1-j2) <= j3
	rule2 = j3 <= j1+j2
	rule3 = m3 == m1 + m2
	rule4 = np.abs(m3) <= j3

	#rules assumed by input
	#assert np.abs(m1) <= j1, 'm1 must be \in {-j1,j1}'
	#assert np.abs(m2) <= j2, 'm2 must be \in {-j2,j2}'

	if rule1 and rule2 and rule3 and rule4:
		#attempting binomial representation
		N1 = np.longdouble((2*j3) + 1 )
		N2 =  np.longdouble(special.factorial(j1 + m1, exact=True)) \
		* np.longdouble(special.factorial(j1 - m1, exact=True)) \
		* np.longdouble(special.factorial(j2 + m2, exact=True)) \
		* np.longdouble(special.factorial(j2 - m2, exact=True)) \
		* np.longdouble(special.factorial(j3 + m3, exact=True)) \
		* np.longdouble(special.factorial(j3 - m3, exact=True))

		N3 = np.longdouble(special.factorial(j1 + j2 - j3, exact=True)) \
		* np.longdouble(special.factorial(j1 - j2 + j3, exact=True)) \
		* np.longdouble(special.factorial(-j1 + j2 + j3, exact=True)) \
		* np.longdouble(special.factorial(j1 + j2 + j3 + 1, exact=True))

		#N = np.longdouble((N1*N2))/np.longdouble((N3))
		N = np.longdouble(0.)
		N += np.divide((N1*N2),N3)

		G = np.longdouble(0.0)

		#k conditions (see eq.5 of https://hal.inria.fr/hal-01851097/document)
		# k  >= 0
		# k <= j1 - m1
		# k <= j2 + m2

		for k in range(0, min([j1-m1, j2+m2]) + 1  ):
			G1 = np.longdouble((-1)**k)
			G2 = np.longdouble(special.comb(j1 + j2 - j3, k,exact=True))
			G3 = np.longdouble(special.comb(j1 - j2 + j3, j1 - m1 - k,exact=True))
			G4 = np.longdouble(special.comb(-j1 +j2 + j3, j2 + m2 - k,exact=True))
			G += np.longdouble(G1*G2*G3*G4)
		Nsqrt = np.longdouble(0)
		Nsqrt += np.sqrt(N)
		return Nsqrt*G

	else:
		return 0.

def clebsch_gordan(l1,m1,l2,m2,l3,m3):
	# try to load c library for calculating cg coefficients
	#if cglib:
	#	return Clebsch_Gordan(l1,m1,l2,m2,l3,m3)
	#else:
	return Clebsch_gordan(l1,m1,l2,m2,l3,m3)

def wigner_3j(j1,m1,j2,m2,j3,m3):
	# uses relation between Clebsch-Gordann coefficients and W-3j symbols to evaluate W-3j
	#VERIFIED - wolframalpha.com
	cg = clebsch_gordan(j1,m1,j2,m2,j3,-m3)

	num = np.longdouble((-1)**(j1-j2-m3))
	denom = np.longdouble(((2*j3) +1)**(1/2))

	return cg* num/denom  # float(num/denom)


def init_clebsch_gordan(lmax):
	#returns dictionary of all cg coefficients to be used at a given value of lmax
	cg = {}
	for l1 in range(lmax+1):
		for l2 in range(lmax+1):
			for l3 in range(lmax+1):
				for m1 in range(-l1,l1+1):
					for m2 in range(-l2,l2+1):
						for m3 in range(-l3,l3+1):
							key = '%d,%d,%d,%d,%d,%d' % (l1,m1,l2,m2,l3,m3)
							cg[key] = clebsch_gordan(l1,m1,l2,m2,l3,m3)
	return cg


def init_wigner_3j(lmax):
	#returns dictionary of all cg coefficients to be used at a given value of lmax
	cg = {}
	for l1 in range(lmax+1):
		for l2 in range(lmax+1):
			for l3 in range(lmax+1):
				for m1 in range(-l1,l1+1):
					for m2 in range(-l2,l2+1):
						for m3 in range(-l3,l3+1):
							key = '%d,%d,%d,%d,%d,%d' % (l1,m1,l2,m2,l3,m3)
							cg[key] = wigner_3j(l1,m1,l2,m2,l3,m3)
	return cg

def store_generalized(coupling_dct,coupling_type,L_R):
	M_Rs = list(coupling_dct.keys())
	ranks = tuple(list(coupling_dct[M_Rs[0]].keys()))
	lmax_per_rank = []
	for rank in ranks:
		max_l = 0
		for lstr in coupling_dct[M_Rs[0]][rank]:
			listr = lstr.split('_')[0]
			lis = [int(k) for k in listr.split(',')]
			maxli = max(lis)
			if maxli > max_l:
				max_l = maxli
		lmax_per_rank.append(max_l)
	lmax_per_rank = tuple(lmax_per_rank)
	ranks_str_lst = ['%d'] *  len(ranks)
	lmax_str_lst = ['%d'] *  len(ranks)
	ranks_str = ''.join(b for b in ranks_str_lst) % tuple(ranks)
	lmax_str = ''.join(b for b in lmax_str_lst) % tuple(lmax_per_rank)
	#coupling_type = 'cg' or 'wig'
	file_name = '%s_LR_%d_r%s_lmax%s.pickle' % (coupling_type,L_R,ranks_str,lmax_str)
	with open(file_name, 'wb') as handle:
		pickle.dump(coupling_dct, handle, protocol=pickle.HIGHEST_PROTOCOL)
	print ('couplings stored in %s' % file_name)

# store a large dictionary of clebsch gordan coefficients
try:
	with open('%s/Clebsch_Gordan.pickle' %lib_path, 'rb') as handle:
		Clebsch_Gordan = pickle.load(handle)
except FileNotFoundError:
	print ("Generating your first pickled library of CG coefficients. This will take a few moments...")
	Clebsch_Gordan = init_clebsch_gordan(lmax_traditional)
	with open('%s/Clebsch_Gordan.pickle' %lib_path, 'wb') as handle:
		pickle.dump(Clebsch_Gordan, handle, protocol=pickle.HIGHEST_PROTOCOL)
# do the same thing for the traditional wigner_3j symbols
try:
	with open('%s/Wigner_3j.pickle' % lib_path, 'rb') as handle:
		Wigner_3j = pickle.load(handle)
except FileNotFoundError:
	print ("Generating your first pickled library of Wigner 3j coefficients. This will take a few moments...")
	Wigner_3j = init_wigner_3j(lmax_traditional)
	with open('%s/Wigner_3j.pickle' % lib_path, 'wb') as handle:
		pickle.dump(Wigner_3j, handle, protocol=pickle.HIGHEST_PROTOCOL)

