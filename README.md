# Permutation-Adapted ACE

To use this module add this directory to your pythonpath:

`
export PYTHONPATH=$PYTHONPATH:/path/to/PA_ACE
`

## Descriptor labels and linear independence

ACE descriptors 


The way the labels for the N-bond ACE descriptors are defined are chemical index of 
the central atom, u0, then chemical indices for other atom in the N bonds, radial 
indices for the N bonds, angular indices for the N bonds, and then intermediate 
angular indices (when applicable). There are N-2 intermediate angular indices for 
descriptors with rank N>=3. For ranks >=3, one will find these intermediate indices 
after the N angular indices for the bonds. In general, the descriptor  labels will 
have the form:
`
u0, u1, u2, ... uN, n1,n2,... nN, l1,l2,... lN, { ... L(N-3), L(N-2) }

For a rank (N=3) example with mu_i=0 , n_i =1, l_i = 0:
0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0
u0, u1, u2, u3, n1, n2, n3, l1, l2, l3, L1

For a non-trivial rank 4 example:
0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0
u0, u1, u2, u3, u4, n1, n2, n3, n4, l1, l2, l3, l4, L1, L2

Where, in the above rank 4 example, only one of two sets of intermediates { (00) , (22) }
yields an independent ACE descriptor. This may be shown analytically using recursion
relationships derived from ladder operations.
`

**The key functionality of this package is to supply users with complete,
linearly independent ACE descriptors analytically as above.**

The labels for linearly independent ACE descriptors may be generated using the following
functions from `pa_gen.py`:


`
from pa_gen import *
pa_labels_raw(rank=4,nmax=4,lmax=4,mumax=1)
`

where `rank` is the number of bonds for the features, `nmax` is the maximum quantum
number for the radial basis, `lmax` is the maximum angular quantum number for the
spherical harmonic basis, and `mumax` is the number of chemical types. This function
will exhaustively generate labels that form a complete, linearly-independent set
of ACE basis functions or descriptors. 


## Coupling coefficients

Both the generalized Wigner symbols and the generalized Clebsch-gordan coefficients
may be generated in this package. The corresponding functions to generate these are given 
in `wigner_couple.py` and `clebsch_couple.py`, respectively and pedagogical implmementations
for calculating them are given in `wigner_tree.py` and `clebsch_tree.py`.

## Pythonized ACE calculator

A stand-alone ACE calculator outside of LAMMPS is provided for testing purposes, 
but is not intended for production usage in simulations. Examples of this may 
be found in `all_examples` and the corresponding functions/classes are given in
`equivariant_calc`.


## LAMMPS usage

This code supplies .yace coupling_coefficient and potential files that may be
ran in LAMMPS. More examples (especially for large lmax and rank) will follow soon.
