from yamlpace_tools.potential import *

# potential parameters. these should match the FitSNAP input
reference_ens = [0.,0.,0.]
# bonds:
bonds=  [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]
#  [(H, H), (H, N), (H, O), (N, N), (N, O), (O, O)]

rcutfaci = [5.0, 5.5, 5.7, 4.4, 5.7, 5.5]
lmbdai = [3.3, 3.3, 3.3, 3.3, 3.3, 3.3]
rcinneri = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
drcinneri = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01]

exhaust_bonds = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)]

rcutfac = []
lmbda = []
rcinner = []
drcinner = []

for bond in exhaust_bonds:
	srt_bnd = tuple(sorted(bond))
	idx = bonds.index(srt_bnd)
	rcutfac.append(rcutfaci[idx])
	lmbda.append(lmbdai[idx])
	rcinner.append(rcinneri[idx])
	drcinner.append(drcinneri[idx])

elements=["H","N","O"]
ranks = [1, 2, 3, 4]
lmax =  [1, 2, 2, 1]
nmax = [2, 1, 1, 1]
lmin = 1
nradbase=max(nmax)


Apot = AcePot(elements,reference_ens,ranks,nmax,lmax,nradbase,rcutfac,lmbda,rcinner,drcinner,lmin,RPI_heuristic='root_SO3_span')
# read the potential file to get expansion coefficients
Apot.write_pot('coupling_coefficients')
