from ..data.monomers import Monomer, Monomers
from ..functions import (
    generate_file,
    PDBtoXYZ,
    CheckMinimalDistance,
)


from unittest import TestLoader, TextTestRunner
from ...tests.unit.test_pdb_file import TestPdbFile
from questionary import text, select, confirm


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
            monomers = Monomers.from_json(monomer_file)
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
            monomers = Monomers.from_json(monomer_file)
            bs_file = text("Enter the name of the monomer's bs file: ").ask()
            monmer = Monomer.create_monomer(bs_file)

            save = confirm("Do you want to save the monomer?").ask()

            monomers.add_monomer(monmer, save=save)
            option = choose_option()

        elif option == "Exit":
            return
