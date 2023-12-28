from raccoon.src.data.monomers import Monomer, Monomers
from raccoon.src.functions import generate_file

monomers = Monomers.from_json('raccoon/src/data/monomers.json')
generate_file(monomers, False, 'gmx/seq_FHFHF.txt', 'gmx/system.pdb')
