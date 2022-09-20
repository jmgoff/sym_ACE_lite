# Example output processing

## FitSNAP generated ACE coefficient file

After FitSNAP fits an ACE potential, a file_name.acecoeff file is written,
containing the expansion coefficients per descriptor. This example shows how
this may be converted to a ctilde potential file that is ready for use in
LAMMPS.

The rewrite_pot.py script calls functions that define the ctilde file,
read in coefficients from the .acecoeff file, and combines them. The inputs
for these functions will reflect inputs for the FitSNAP run. The FitSNAP 
infile is provided for reference as well as a summary for an example fit with UO2.


