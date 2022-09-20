from rpi_lib import *
import sys

def get_munl_tups(nus):
	unique_munls = []
	for nu in nus:
		mu0,mu,n,l,L = get_mu_n_l(nu,return_L=True)
		munl = tuple([tuple(mu),tuple(n),tuple(l)])
		if munl not in unique_munls:
			unique_munls.append(munl)
	return unique_munls

def arrange_by_degree(nus,return_degree=False,print_info=False):
	degrees = []
	deg_per_nu = {nu:None for nu in nus}
	for nu in nus:
		mu0,mu,n,l,L = get_mu_n_l(nu,return_L=True)
		deg = np.sum(mu) + np.sum(n) + np.sum([i for i in l])
		deg_per_nu[nu] = deg
		degrees.append(deg)
	unique_degs = sorted(list(set(degrees)))

	nu_per_deg = {deg:[] for deg in unique_degs}
	for nu in nus:
		deg = deg_per_nu[nu]
		nu_per_deg[deg].append(nu)

	if print_info:
		print (nu_per_deg)

	deg_nu_sublists = [ nu_per_deg[deg] for deg in unique_degs]
	if return_degree:
		return deg_nu_sublists, unique_degs
	else:
		return deg_nu_sublists

rank = int(sys.argv[1])
nmax = int(sys.argv[2])
lmax = int(sys.argv[3])
mumax = int(sys.argv[4])
try:
	lmin = int(sys.argv[5])
except IndexError:
	lmin=0

import time
t1 = time.time()
symmetric_set = descriptor_labels_YSG(rank,nmax,lmax,mumax,lmin)
t2 = time.time()
# flag to print the counts for lexicographically ordered starting set for SVD
print_lexico=True
#flag to print polynomial degree (True) or scaled degree (False)
print_pure=True
#flag to print munlL labels for input to lammps
print_all_labels=False
#flag to print cumulative vs exact descriptor counts per degree
cumulative = True
if print_lexico:
	lexicographic_set = generate_nl(rank,nmax,lmax,mumax=mumax,lmin=lmin,all_perms=True)
	lexicographic_set2 = generate_nl(rank,nmax,lmax,mumax=mumax,lmin=lmin,all_perms=False)
sym_deg,degrees = arrange_by_degree(symmetric_set,return_degree=True)
if print_lexico:
	lex_deg = arrange_by_degree(lexicographic_set)
	lex_deg2 = arrange_by_degree(lexicographic_set2)
if print_all_labels:
	for nl in symmetric_set:
		print ('"%s",' % nl)

sym_deg_lens = [len(b) for b in sym_deg]
if print_lexico:
	lex_deg_lens= [len(b) for b in lex_deg]
	lex_deg_lens2= [len(b) for b in lex_deg2]

if cumulative:
	max_deg_per_rank = (rank/max(degrees))
	
	if print_pure:
		print ( ' , '.join(str(b) for b in degrees))
		print ( ' , '.join(str(b) for b in np.cumsum(sym_deg_lens)))
		if print_lexico:
			print ( ' , '.join(str(b) for b in np.cumsum(lex_deg_lens2)))
			print ( ' , '.join(str(b) for b in np.cumsum(lex_deg_lens)))
	else:
		print ( ' , '.join(str(b/max_deg_per_rank) for b in degrees))
		print ( ' , '.join(str(b) for b in np.cumsum(sym_deg_lens)))
		if print_lexico:
			print ( ' , '.join(str(b) for b in np.cumsum(lex_deg_lens2)))
			print ( ' , '.join(str(b) for b in np.cumsum(lex_deg_lens)))
else:
	print ( ' '.join(str(b) for b in degrees))
	print ( ' '.join(str(b) for b in sym_deg_lens))
	if print_lexico:
		print ( ' '.join(str(b) for b in lex_deg_lens))

print (t2-t1)
