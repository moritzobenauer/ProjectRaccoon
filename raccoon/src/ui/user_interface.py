from ..data.monomers import Monomer, Monomers
from ..functions import (
    generate_file,
    get_atoms_from_bs_file,
    PDBtoXYZ,
    CheckMinimalDistance,
)

from unittest import TestLoader, TextTestRunner, TestCase
from ...tests.unit.test_pdb_file import TestPdbFile
from questionary import text, select, confirm


def add_monomer(fpath: str):
    """Adds a monomer to the monomers.dat file.

    Args:
        fpath (str): input file
    """

    # create empty monomer
    monomer = Monomer()

    atoms = get_atoms_from_bs_file(fpath)

    for idx, atom in enumerate(atoms):
        identifier = text(
            f"Enter a force field Identifier for {atom[1]}' with the Number {idx}: "
        ).ask()
        atoms.append(identifier)

    monomer.atoms = atoms

    monomer.name = text("Enter the name of the monomer").ask()
    monomer.polymer = confirm("Is this a standard polymer?").ask()
    monomer.resolution = select(
        "Choose resolution",
        choices=["atomistic", "united_atom", "coarse_grained"],
    ).ask()

    linkC = text(f"Choose C-Terminus (1-{len(atoms)})").ask()
    linkN = text(f"Choose N-Terminus (1-{len(atoms)})").ask()
    monomer.link = [int(linkC), int(linkN)]

    monomer.inverted = False

    if confirm("Save the monomer to the monomer data file?").ask():
        monomer.add_to_file()

    return


def choose_option() -> str:
    return select(
        "Choose Function",
        choices=[
            "Create PDB File",
            "Check PDB File",
            "Convert PDB to XYZ File",
            "Check Minimal Distance",
            "Exit",
        ],
    ).ask()


def start_racoon(
    sequence_file: str,
    out_file: str,
    monomer_file: str,
    explicitbonds: bool,
    remove_duplicates: bool,
):
    option = choose_option()

    while True:
        if option == "Create PDB File":
            monomers = Monomers.from_file(monomer_file)
            generate_file(monomers, explicitbonds, sequence_file, out_file)

            option = choose_option()

        elif option == "Check PDB File":
            suite = TestLoader().loadTestsFromTestCase(TestPdbFile)
            TextTestRunner().run(suite)
            option = choose_option()

        elif option == "Convert PDB to XYZ File":
            PDBtoXYZ(out_file)
            option = choose_option()

        elif option == "Check Minimal Distance":
            CheckMinimalDistance(out_file)
            option = choose_option()

        elif option == "Add Monomer":
            inputfile = text("Enter the path to the monomer file").ask()
            add_monomer(inputfile)
            option = choose_option()

        elif option == "Exit":
            return
