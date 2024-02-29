from project_raccoon.src.data.monomers import Monomer, Monomers
from project_raccoon.src.functions import generate_file

monomers = Monomers.from_file("monomers.dat")

generate_file(monomers, False, "gmx/seq_FHFHF.txt", "gmx/system.pdb", trr=1)
