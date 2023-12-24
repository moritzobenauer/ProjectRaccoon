from raccoon.src.data.monomers import Monomer, Monomers

from collections import namedtuple
from raccoon.src.typing import List, Dict, Tuple, NamedTuple

import numpy as np


def generate_sequence(
    monomers: Monomers, fpath: str
) -> NamedTuple("sequence", [("index", List), ("inverted", List), ("reps", List)]):
    """
    Generates a sequence from a inputfile.

    Args:
        monomers (Monomers): Monomers object.
        fpath (str): Path to the file.

    Returns:
        namedtuple: Sequence which contains the index of the monomer, the information if it is inverted and the number of repetitions.
    """

    sequence = namedtuple("sequence", ["index", "inverted", "reps"])
    sequence.index = list()
    sequence.inverted = list()
    sequence.reps = list()

    with open(fpath, "r") as f:
        for line in f.readlines():
            if line.startswith("#") or line.strip() == "":
                continue
            else:
                res, resolution, inverted, reps = line.split(":")
                if resolution == "AA":
                    resolution_lookup = "atomistic"
                elif resolution == "UA":
                    resolution_lookup = "united_atom"
                elif resolution == "CG":
                    resolution_lookup = "coarse_grained"

                index = monomers.index({"name": res, "resolution": resolution_lookup})

                sequence.index.append(int(index))
                sequence.inverted.append(bool(inverted))
                sequence.reps.append(int(reps))

    return sequence


def generate_file(monomers: Monomers, explicit_bonds: bool, spath: str, outpath: str):
    """Central function of the modul: adds monomers to a polymer peptide chain and writes it to a PDB file.

    Args:
        monomers (Monomers):
        explicit_bonds (bool):
        outpath (str):
        spath (str):
    """
    atom_count = 0
    res_count = 0

    with open(outpath, "w") as f:
        x_min, x_max = -5.5, 5.5
        y_min, y_max = -5.5, 5.5
        z_min, z_max = 0.3, 0.4

        # cartesian shifts
        cshifts = np.zeros(3)
        atom_count = 0

        links = list()

        sequence = generate_sequence(monomers, spath)
        for index, inverted, reps in zip(
            sequence.index, sequence.inverted, sequence.reps
        ):
            monomer = monomers[index]

            if inverted:
                monomer = monomer.invert()

            # TODO: implement explicit bonds
            # if explicit_bonds:
            #    ExplicitBonds(monomers[monomer.index])

            for rep in range(reps):
                m = np.array(
                    [
                        (np.random.uniform(x_min, x_max)),
                        (np.random.uniform(y_min, y_max)),
                        (float(monomer.atom_count) * np.random.uniform(z_min, z_max)),
                    ]
                )

                cshifts += m  # TODO: += oder = m? soll der shift mit der anzahl der atome größer werden?

                updated_monomer = monomer.update(atom_count, cshifts)

                atom_count += updated_monomer.atom_count
                res_count += 1

                # add the C-Terminus link and / or the N-Terminus link
                links.extend(updated_monomer.link)

                for atom in updated_monomer.atoms:
                    f.write(
                        "{:>0}{:<7}{:<5}{:<5}{:<4}{:<3}{:<6}{:<8}{:<8}{:<10}{:<7}{:<14}{}\n".format(
                            "",
                            "ATOM",
                            atom[6],
                            atom[0],
                            updated_monomer.name,
                            "A",
                            res_count,
                            atom[2],
                            atom[3],
                            atom[4],
                            1.0,
                            0.0,
                            atom[1],
                        )
                    )

                # TODO: WriteExplicitBonds(out_file)

        # write bonds in file
        for i in range(0, len(links) - 1, 1):
            f.write(
                "{:>0}{:<7}{:<5}{:<5}{:<4}{:<3}{:<6}{:<8}{:<8}{:<10}{:<7}{:<14}{}\n".format(
                    "",
                    "CONECT",
                    links[i],
                    links[i + 1],
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                )
            )

    # no need for sorting beacuse the atoms are already sorted ( consecutive numbering) and the bonds are written after all atoms
    # sort_PDB(outpath)
    close_PDB(outpath, atom_count)

    print(
        f"Created PDB file with {res_count} residue and {atom_count} atoms to {outpath}."
    )


def sort_PDB(fpath: str):
    """
    Sorts a PDB file by the atom index.

    Args:
        file_path (str): Path to the file.

    """

    with open(fpath, "r+") as f:
        rows = f.readlines()
        atoms = [row for row in rows if row.startswith("ATOM")]
        bonds = [row for row in rows if row.startswith("CONECT")]

        sorted_atoms = sorted(atoms, key=lambda x: x[1])
        sorted_bonds = sorted(bonds, key=lambda x: x[1])

    with open(fpath, "w") as f:
        f.writelines(sorted_atoms)
        f.writelines(sorted_bonds)


def close_PDB(fpath: str, atom_count: int):
    """
    Closes a PDB file by adding the END statement and the number of atoms.

    Args:
        file_path (str): Path to the file.
        atom_count (int): Number of atoms.

    """

    with open(fpath, "a") as file:
        file.write(
            "{:>0}{:<11}{:<5}{:<5}{:<4}{:<5}{:<5}{:<5}{:<5}{:<5}{:<5}{:<5}{:<5}{}\n".format(
                "", "MASTER", 0, 0, 0, 0, 0, 0, 0, 0, atom_count, 0, atom_count, 0
            )
        )
        file.write("END")
