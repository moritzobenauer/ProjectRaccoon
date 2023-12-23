from ..data.monomers import Monomers
from ..functions import generate_file


def start_racoon(
    sequenze_file: str,
    monomer_file: str,
    out_file: str,
    explicitbonds: bool,
    remove_duplicates: bool,
):
    global monomers

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
        elif options == "Exit":
            return
