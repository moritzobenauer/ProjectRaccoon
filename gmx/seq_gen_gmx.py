from raccoon.src.data.monomers import Monomer, Monomers
from raccoon.src.functions import generate_file

monomers = Monomers.from_file('monomers.dat')


for k in range(1000):
    # Polypeptide sequences with up to 350 atoms / beads work realiably with trr <= 0.5
    for t in range(1,5,1):

        q = t/10
        # Polypeptides FHFHF Sequence
        generate_file(monomers, False, 'gmx/seq_FHFHF.txt', 'gmx/system.pdb', trr=q)
         # Polymer (PEO)50
        generate_file(monomers, False, 'gmx/seq_PEO50.txt', 'gmx/peo50.pdb', trr=q)
        # Polymer Peptide Conjugate
        generate_file(monomers, False, 'gmx/seq_AcFGPEO50GFAc.txt', 'seq_AcFGPEO50GFAc.pdb', trr=q)

    # Linear polymer sequences with up to 700 atoms / beads work reliably with trr <= 0.6

    for t in range(1,6,1):
        q = t/10
        # Polymer (PEO)100
        generate_file(monomers, False, 'gmx/seq_PEO100.txt', 'gmx/peo100.pdb', trr=q)

generate_file(monomers, False, 'gmx/seq_FHFHF.txt', 'gmx/system.pdb', trr=1)