from raccoon.src.data import Monomer, Monomers, Sequence, Atom

from collections import namedtuple
from raccoon.src.typing import List, Dict, Tuple, NamedTuple

import numpy as np


def generate_sequence(monomers: Monomers, fpath: str) -> Sequence:
    """
    Generates a sequence from a inputfile.

    Args:
        monomers (Monomers): Monomers object.
        fpath (str): Path to the file.

    Returns:
        namedtuple: Sequence which contains the index of the monomer, the information if it is inverted and the number of repetitions.
    """

    index = list()
    inverted = list()
    reps = list()

    with open(fpath, "r") as f:
        for line in f.readlines():
            if line.startswith("#") or line.strip() == "":
                continue
            else:
                res, resolution, inv, rep = line.split(":")
                if resolution == "AA":
                    resolution_lookup = "atomistic"
                elif resolution == "UA":
                    resolution_lookup = "united_atom"
                elif resolution == "CG":
                    resolution_lookup = "coarse_grained"

                index.append(
                    monomers.index({"name": res, "resolution": resolution_lookup})
                )
                inverted.append(bool(inv))
                reps.append(int(rep))

    return Sequence(index, inverted, reps)


def RandShift(
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    z_min: float,
    z_max: float,
    z_bias: float,
):
    """Creates a random shift in all directions with and without bias. Setting r_(min, _max) = 1 creates a linear shift for debugging. No default values needed.

    Args:
        x_min (float): _description_
        x_max (float): _description_
        y_min (float): _description_
        y_max (float): _description_
        z_min (float): _description_
        z_max (float): _description_
        z_bias (float): _description_

    Returns:
        array (np.array): 3d vector with shape (3,)
    """
    x = np.random.uniform(x_min, x_max)
    y = np.random.uniform(y_min, y_max)
    z = z_bias * np.random.uniform(z_min, z_max)
    return np.array([x, y, z])


def MinimalDistance(
    polypetide_coordinates: np.array, new_monomer_coordinates: np.array
):
    """Checks for the minimal distance between atoms / beads of the new monomer and all prior atoms / beads.
       Future revisions could use contacts matrix to adjust self-avoiding random walk.

    Args:
        coordinates (np.array): Array of coordinates of the prior atoms / beads.
        new_monomer (np.array): Array of coordinates of the new monomer.

    Returns:
        min (float): Minimal distance between atoms / beads.
        contacts (np.array): Matrix of distances between atoms / beads.
    """
    contacts = np.empty(
        (polypetide_coordinates.shape[0], new_monomer_coordinates.shape[0])
    )
    for i, r1 in enumerate(new_monomer_coordinates):
        for j, r2 in enumerate(polypetide_coordinates):
            contacts[j, i] = np.round(np.linalg.norm(r1 - r2), 3)
    min = np.min(contacts)
    return min, contacts


# Self-Avoiding Random Walk to prevent infinite forces upon energy minimization
# Combines RandShift() and MinimalDistance(). Returns k, which can be used to add onto the new monomer.
# Treshshold of trr=1 should be sufficient to prevent infinite forces


def SemiRandomWalk(
    polypeptide_coordinates: np.array, monomer: Monomer, trr: float, shift: List[float]
):
    """Self-Avoiding Random Walk to prevent infinite forces upon energy minimization.
         Combines RandShift() and MinimalDistance(). Returns k, which can be used to add onto the new monomer.
            Treshshold of trr=1 should be sufficient to prevent infinite forces.

    Args:
        polypetide_coordinates (np.array): Array of all coordinates of the atoms in the previous monomers.
        monomer (Monomer): Monomer from raccoon.src.data
        trr (float): Threshold for minimal distance.
        shift (list(float)): List of cartesian shifts.

    Returns:
        k (np.array): 3d vector with shape (3,)
    """

    minimal_distance = 0
    monomer_coordinates = monomer.coordinates_to_numpy()
    while minimal_distance < trr:
        k = RandShift(*shift)
        updated_monomer_coordinates = monomer_coordinates + k[np.newaxis, :]
        minimal_distance, _ = MinimalDistance(
            polypeptide_coordinates, updated_monomer_coordinates
        )
    return (updated_monomer_coordinates - monomer_coordinates)[0]


def generate_file(
    monomers: Monomers,
    explicit_bonds: bool,
    spath: str,
    outpath: str,
    trr: float = 1,
    shift_cartesian: List[float] = [-1, 1, -1, 1, -1, 1, 1],
    damping_factor: float = 0.5,
):
    """Central function of the modul: adds monomers to a polymer peptide chain and writes it to a PDB file.

    Args:
        monomers (Monomers):
        explicit_bonds (bool):
        outpath (str):
        spath (str):
        trr (float):
        shift_cartesian (list(float)):
        damping_factor (float):
    """
    atom_count = 0
    res_count = 0

    with open(outpath, "w") as f:
        # cartesian shifts
        cshifts = np.zeros(3)
        atom_count = 0

        links = list()
        links_explicit = list()
        coordinates = np.zeros((1, 3))

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
                shift_cartesian[6] = float(monomer.atom_count) * damping_factor
                m = SemiRandomWalk(coordinates, monomer, trr=1, shift=shift_cartesian)

                cshifts += m

                updated_monomer = monomer.update(atom_count, cshifts)

                new_row = monomer.coordinates_to_numpy()
                coordinates = np.vstack((coordinates, new_row))

                pairs = []

                # Probably not the most efficient method to get all the explicit links, it works however
                # updated_monomer.get_explicit_links() gets a list of all links of every atom in a monomer
                # iterate through that list and create pairs

                if explicit_bonds == True:
                    for index, neighbor in enumerate(
                        updated_monomer.get_explicit_links()
                    ):
                        for n in neighbor:
                            pairs.append((index + 1 + atom_count, n))
                    unique_pairs = set()

                    # search for duplicate touples and only keep the non-duplicates

                    for item in pairs:
                        if (item not in unique_pairs) and (
                            tuple(reversed(item)) not in unique_pairs
                        ):
                            unique_pairs.add(item)

                    for items in unique_pairs:
                        links_explicit.extend(items)

                else:
                    pass

                atom_count += updated_monomer.atom_count
                res_count += 1

                # add the C-Terminus link and/or the N-Terminus link
                links.extend(updated_monomer.link)

                for atom in updated_monomer.atoms:
                    f.write(
                        "{:>0}{:<7}{:<5}{:<5}{:<4}{:<3}{:<6}{:<8}{:<8}{:<10}{:<7}{:<14}{}\n".format(
                            "",
                            "ATOM",
                            atom.index,
                            atom.ff_identifier,
                            updated_monomer.name,
                            "A",
                            res_count,
                            atom.x,
                            atom.y,
                            atom.z,
                            1.0,
                            0.0,
                            atom.element,
                        )
                    )

        # extra loop needed, since iterator has to be set to [::2] instead of [::1] for non-explicit bonds

        if explicit_bonds == True:
            for i in range(0, len(links_explicit) - 1, 2):
                f.write(
                    "{:>0}{:<7}{:<5}{:<5}{:<4}{:<3}{:<6}{:<8}{:<8}{:<10}{:<7}{:<14}{}\n".format(
                        "",
                        "CONECT",
                        links_explicit[i],
                        links_explicit[i + 1],
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
