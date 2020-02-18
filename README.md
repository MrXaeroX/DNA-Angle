# DNA-Angle
DNA bending angle calculator for arbitrary PDB structures.
Automatically detects complimentary residues, builds control points for
the phosphate backbone, approximates them with a polynomial regression,
and outputs angle (in degrees) between the average tangents at the ends
of the DNA. Requires NumPy to run.