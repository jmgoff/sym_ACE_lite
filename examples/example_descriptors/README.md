# Examples

## Descriptor calculations
ACE descriptors may be calculated for testing purposes using the libraries in 
descriptor_calc_local. The examples below show how a set of descriptors may be
calculated for a given atomic structure. 

For production MD, fitting of potentials, or processing of large datasets, it
is recommended to use the LAMMPS compute.

## Usage:
<code>
from descriptor_calc_local.chem_B import *  
  
ranks = range(1 , 6 + 1)                                # ranks of basis functions to be evaluated<br/>
elements = ['Ag']                                       # list of elements<br/>
rc = 20.5                                               # radial cutoff (angstroms) for radial basis<br/>
lmbda = 5.0                                             # exponential factor for sampling of r << r_c in g_k(r)<br/>
nradbase = 16                                           # For now, this is just the n_max for rank = 1<br/>
lmax_dict = {1:0,2:2,3:2,4:2,5:2,6:2}                   # maximum l quantum number in basis functions per rank: {rank:lmax(rank)} (note lmax of rank 1 must be 0 for symmetry reduction)<br/>
nradmax_dict = {1:nradbase, 2: 2, 3:2, 4:2, 5:2, 6:2}   # maximum n quantum number in basis functions per rank: {rank:nmax(rank)}<br/>
i_d = 0                                                 # index supplied to identify the atomic configuration (useful when iterating over many configurations)<br/>
</code>

Parameters are specified as in previous versions of Wigner_ACE. Then descriptor labels 
are generated using the lexicographical ordering scheme (Drautz 2019, 2020 or Dusson 2022)
or the symmetry-reduced scheme. The functions are 'generate_nl' and 'descriptor_labels', 
respectively to generate the corresponding sets of descriptor labels. It is noted that
no SVD is performed on the lexicographically ordered set as is done in Dusson 2022 or 
Drautz 2020.

An ASE atoms object is then constructed, and descriptors may be calculated using internal
python code using the 'descriptor_calc_local.chem_b' module. The 'atoms_eB' function in this 
module computes a list of descriptors for a given ASE atoms object. The output is a dictionary
over each atom in the structure, providing all descriptor values per atom.

<code>
```
\# input for descriptor calculator is expected in dictionary format  
arg = {'id':i_d,'atoms':atoms,'elements':elements,'rc':rc,'nradbase':nradbase,  
'nradmax':nradmax_dict,'lmax':lmax_dict,'lmbda':lmbda,'rank':ranks,'nus':nus, 'L_R':0, 'M_R':0}  

B_per_atom = atoms_eB(arg)  

import json  
with open('outfile.json','w') as writejson:  
        json.dump(B_per_atom, writejson, sort_keys=False, indent=2)  
```
</code>

For the set of descriptors coming from lexicographical ordering, the example 
'Ag_lexico_base.py' is provided. Running this will produce an output file containing all
descriptor values atoms in a random silver cluster called 'lexico_ag.json'. The example
'Ag_symmetric_base.py', does the same, but using the symmetry reduced descriptor set.
The corresponding file is 'symmetric_ag.json'. Example outputs are provided in 
\<path_to\>/symmetric_ACE/examples/ref, and should be reproducable up to machine precision.
