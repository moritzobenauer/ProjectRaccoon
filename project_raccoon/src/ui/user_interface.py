from ..data.monomers import Monomer, Monomers
from ..functions import (
    generate_sequence,
    generate_file,
    pdb_to_xyz,
    check_minimal_distance,
    check_pdb_file,
)

from rich.console import Console
from questionary import text, select, confirm


def choose_option() -> str:
    return select(
        "Choose Function",
        choices=[
            "Create PDB File",
            "Check PDB File",
            "Convert PDB to XYZ File",
            "Check Minimal Distance",
            "Manage Monomers",
            "Exit",
        ],
    ).ask()


def manage_monomers() -> str:
    return select(
        "Choose Function",
        choices=[
            "Add Monomer",
            "Delete Monomer",
            "Print Monomers",
            "Export JSON Monomer File",
            "Return",
        ],
    ).ask()


def start_racoon(
    sequence_file: str,
    out_file: str,
    monomer_file: str,
    explicitbonds: bool,
    remove_duplicates: bool,
    suppress_messages: bool = False,
):
    option = choose_option()

    monomers = Monomers.from_json(monomer_file)

    # Use rich terminal to make it look a bit nicer :)
    console = Console()

    try:
        while True:
            if option == "Create PDB File":
                try:
                    sequence = generate_sequence(monomers=monomers, fpath=sequence_file)
                    generate_file(
                        monomers=monomers,
                        sequence=sequence,
                        explicit_bonds=explicitbonds,
                        outpath=out_file,
                        suppress_messages=suppress_messages,
                    )
                except Exception as e:
                    import sys

                    print(
                        "Caught the following error while generating the PDB file:",
                        file=sys.stderr,
                    )
                    print(
                        "---------------------------------------------------------",
                        file=sys.stderr,
                    )
                    print(e, file=sys.stderr)
                    print(
                        "---------------------------------------------------------",
                        file=sys.stderr,
                    )
                    print("Exiting", file=sys.stderr)
                    exit(1)

                option = choose_option()

            elif option == "Check PDB File":
                try:
                    check_pdb_file(out_file, suppress_messages=suppress_messages)
                except Exception as e:
                    console.print(f"Error: {e}", style="bold red")

                option = choose_option()

            elif option == "Convert PDB to XYZ File":
                pdb_to_xyz(out_file, suppress_messages=suppress_messages)
                option = choose_option()

            elif option == "Check Minimal Distance":
                check_minimal_distance(out_file)
                option = choose_option()

            elif option == "Manage Monomers":
                sec_option = manage_monomers()

                if sec_option == "Print Monomers":
                    console.print(
                        [
                            f"{index} {monomer.name} {monomer.resolution}"
                            for index, monomer in enumerate(monomers)
                        ]
                    )
                    sec_option = manage_monomers()

                elif sec_option == "Delete Monomer":
                    string_monomers = [
                        [
                            f"{index} {monomer.name} {monomer.resolution}"
                            for index, monomer in enumerate(monomers)
                        ]
                    ][0]
                    monomer_select = select(
                        "Choose Monomer to remove", choices=string_monomers
                    ).ask()
                    monomer_identifier = int(monomer_select.strip()[0])

                    save = confirm(
                        f"Do you want to remove the monomer {monomer_select}?"
                    ).ask()
                    if save:
                        monomers.remove_monomer(
                            monomer=monomers[monomer_identifier], save=save
                        )
                        console.print(
                            f"{monomer_select} was removed.", style="bold red"
                        )
                    if not save:
                        console.print("No monomers were removed.", style="bold red")

                elif sec_option == "Export JSON Monomer File":
                    fpath = text("Enter JSON output filename").ask()
                    if not fpath.split(".")[-1] == "json":
                        console.print(
                            "Please specify a JSON output file.", style="bold red"
                        )
                        return
                    monomers.to_json(fpath)
                    sec_option = manage_monomers()

                elif sec_option == "Add Monomer":
                    name = text("Enter the name of the monomer").ask()
                    while name == "":
                        console.print("Please enter a valid name.", style="bold red")
                        name = text("Enter the name of the monomer").ask()
                    resolution = select(
                        "Choose resolution",
                        choices=["atomistic", "united_atom", "coarse_grained"],
                    ).ask()
                    polymer = confirm("Is this a polymer?").ask()

                    bs_file = text("Enter the name of the monomer's bs file: ").ask()
                    try:
                        atoms = Monomer.get_atoms_from_bs_file(bs_file)
                    except Exception as e:
                        console.print(f"Error: {e}", style="bold red")
                        sec_option = manage_monomers()

                    [
                        console.print(
                            f"Index: {index+1} Element Symbol: {atom[0]} Neighbors: {atom[4]}"
                        )
                        for index, atom in enumerate(atoms)
                    ]

                    linkC = text(f"Choose C-Terminus (1-{len(atoms)})").ask()
                    linkN = text(f"Choose N-Terminus (1-{len(atoms)})").ask()
                    link = [int(linkC), int(linkN)]

                    ff_identifiers = list()
                    for atom in atoms:
                        ff_identifier = text(
                            f"Enter a force field Identifier for {atom[0]}' with the Number {atom[-1]}: "
                        ).ask()
                        ff_identifiers.append(ff_identifier)

                    monmer = Monomer.create_monomer(
                        name,
                        resolution,
                        polymer,
                        link,
                        atoms,
                        ff_identifiers,
                    )

                    save = confirm("Do you want to save the monomer?").ask()

                    monomers.add_monomer(monmer, save=save)

                    sec_option = manage_monomers()

                elif sec_option == "Return":
                    option = choose_option()

            elif option == "Exit":
                if not suppress_messages:
                    console.print("Project RACCOON was successfully finished.")
                return

    except KeyboardInterrupt:
        console.print("Project RACCOON cancelled by user.", style="bold red")
        return
