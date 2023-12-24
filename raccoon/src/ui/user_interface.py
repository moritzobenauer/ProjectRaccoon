from ..data.monomers import Monomer, Monomers
from ..functions import generate_file, get_atoms_from_bs_file

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


def start_racoon(
    sequenze_file: str,
    monomer_file: str,
    out_file: str,
    explicitbonds: bool,
    remove_duplicates: bool,
):
    options = select(
        "Choose Function",
        choices=[
            "Create PDB File",
            "Check PDB File",
            "Convert PDB to XYZ File",
            "Check Minimal Distance",
            "Exit",
        ],
    ).ask()

    while True:
        if options == "Create PDB File":
            monomers = Monomers.from_file(monomer_file)
            generate_file(monomers, explicitbonds, sequenze_file, out_file)

        # elif options == 'Check PDB File':
        #    #CheckPDB(outputfile)
        #    #main()
        # elif options == "Convert PDB to XYZ File":
        #    #PDBtoXYZ(outputfile)
        #    #main()
        # elif options == 'Check Minimal Distance':
        #    #CheckMinimalDistance(outputfile)
        #    #main()
        elif options == "Add Monomer":
            # create empty monomer

            inputfile = text("Enter the path to the monomer file").ask()

            add_monomer(inputfile)

        elif options == "Exit":
            return
