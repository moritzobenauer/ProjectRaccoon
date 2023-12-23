from raccoon.src.data.monomers import Monomer, Monomers

from collections import namedtuple
from typing import List, Dict, Tuple

import numpy as np


def generate_sequence(monomers: Monomers, fpath: str) -> List[int]:
    """
    Generates a sequence from a inputfile.

    Args:
        monomers (Monomers): Monomers object.
        fpath (str): Path to the file.

    Returns:
        namedtuple: Sequence which contains the index of the monomer, the information if it is inverted and the number of repetitions.
    """

    sequence = namedtuple("sequence", ["index", "inverted", "repeat"])
    sequence.index = list()
    sequence.inverted = list()
    sequence.reps = list()

    with open(fpath, "r") as input:
        for line in input:
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

                index = monomers.index({"res": res, "resolution": resolution_lookup})

                sequence.index.append(index)
                sequence.inverted.append(inverted)
                sequence.reps.append(reps)

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
        ashift = 0

        links = list()

        for index, inverted, reps in generate_sequence(monomers, spath):
            monomer = monomers[index]
            if inverted == "True":
                monomer = monomer.invert()

            # TODO: implement explicit bonds
            # if explicit_bonds:
            #    ExplicitBonds(monomers[monomer.index])

            for rep in range(reps):
                m = np.array(
                    [
                        (np.random.uniform(x_min, x_max)),
                        (np.random.uniform(y_min, y_max)),
                        (monomer.atom_count * np.random.uniform(z_min, z_max)),
                    ]
                )

                cshifts += m  # TODO: += oder = m? soll der shift mit der anzahl der atome größer werden?

                updated_monomer = monomer.update(atom_count, cshifts)

                atom_count += monomer.atom_count
                res_count += 1

                # add the C-Terminus link
                links.append(updated_monomer.link[1])
                # add the N-Terminus link
                links.append(updated_monomer.link[0])

                for atom in updated_monomer.atoms:
                    f.write(
                        "{:>0}{:<7}{:<5}{:<5}{:<4}{:<3}{:<6}{:<8}{:<8}{:<10}{:<7}{:<14}{}\n".format(
                            "",
                            "ATOM",
                            atom[6],
                            atom[0],
                            updated_monomer["res"],
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
                for i in range(0, len(links), 1):
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
    sort_PDB(outpath)
    close_PDB(outpath, atom_count)

    print(f"Created PDB file with {res_count} residue and {atom_count} atoms.")


def sort_PDB(fpath: str):
    """
    Sorts a PDB file by the atom index.

    Args:
        file_path (str): Path to the file.

    """

    with open(fpath, "r") as f:
        rows = f.readlines()
        sorted_rows = sorted(rows, key=lambda x: x.split()[0])
        f.writelines(sorted_rows)


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
