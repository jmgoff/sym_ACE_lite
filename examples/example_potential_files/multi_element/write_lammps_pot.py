from yamlpace_tools.potential import *


elements=['In','P']
reference_ens=[-1.65967588701534, 4.38159549501534]
nmax=[16,2,2,2]
lmax=[0,2,2,2]
ranks=[1,2,3,4]
nradbase=max(nmax)
rcut = {'[0, 1]': 7.5, '[1, 0]': 7.5, '[1, 1]': 5.0, '[0, 0]': 4.5}
lmbda = {'[0, 1]': 3.0, '[1, 0]': 3.0, '[1, 1]': 2.0, '[0, 0]': 2.0}
rcutinner = {'[0, 1]': 0.5, '[1, 0]': 0.5, '[1, 1]': 0.6, '[0, 0]': 0.4}
drcutinner = {'[0, 1]': 0.01, '[1, 0]': 0.01, '[1, 1]': 0.01, '[0, 0]': 0.01}

#permutation symmetry adapted ACE labels
Apot = AcePot(elements,reference_ens,ranks,nmax,lmax,nradbase,rcut,lmbda,rcutinner,drcutinner,RPI_heuristic='root_SO3_span')

#Full lexicographical set (no numerical reduction)
#Apot = AcePot(elements,reference_ens,ranks,nmax,lmax,nradbase,rcut,lmbda,rcutinner,drcutinner,RPI_heuristic='lexicographical')

#Write coupling coefficient file
#This is the file needed to run the lammps compute
Apot.write_pot('coupling_coefficients')

# If expansion coefficients are known, they may be added to the file and a potential may be written
# NOTE that elements should be sorted alphabetically for your coefficients
#betas = np.random.rand(len(nus))
#Apot.set_betas(betas)
# Reset functions with added expansion coefficients - note that a coupling coefficient file may not be written after expansion coefficients are set.
#Apot.set_funcs()
#Apot.write_pot('potential')
