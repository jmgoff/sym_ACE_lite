from yamlpace_tools.potential import *

# Write ML-PACE ready potential files
elements=['Mo']
reference_ens=[0.]
nmax=[16,2,2,2]
lmax=[0,2,2,2]
ranks=[1,2,3,4]
nradbase=max(nmax)
rcut=7.5
lmbda=3.33

#permutation symmetry adapted ACE labels
Apot = AcePot(elements,reference_ens,ranks,nmax,lmax,nradbase,rcut,lmbda,RPI_heuristic='root_SO3_span')

#Full lexicographical set (no numerical reduction)
#Apot = AcePot(elements,reference_ens,ranks,nmax,lmax,nradbase,rcut,lmbda,RPI_heuristic='lexicographical')

#Write coupling coefficient file
#This is the file needed to run the lammps compute
Apot.write_pot('coupling_coefficients')

# If expansion coefficients are known, they may be added to the file and a potential may be written
#betas = np.random.rand(len(nus))
#Apot.set_betas(betas)
# Reset functions with added expansion coefficients - note that a coupling coefficient file may not be written after expansion coefficients are set.
#Apot.set_funcs()
#Apot.write_pot('potential')
