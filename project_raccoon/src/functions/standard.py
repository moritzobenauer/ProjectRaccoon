from ..data import Monomer, Monomers, Sequence, Atom

from ..typing import List, Dict, Optional
from ..util import eps

import numpy as np
from rich.console import Console
from tqdm import tqdm


from scipy.spatial.distance import cdist


def generate_sequence(monomers: Monomers, fpath: str) -> Sequence:
    """
    Generates a sequence from a inputfile.

    Args:
        - `monomers (Monomers)`: Monomers object.
        - `fpath (str)`: Path to the file.

    Returns:
        - `sequence (Sequence)`: Sequence which contains the index of the monomer, the information if it is inverted and the number of repetitions.
    """

    index = list()
    inverted = list()
    reps = list()

    with open(fpath, "r") as f:
        for i, line in enumerate(f.readlines()):
            if line.startswith("#") or line.strip() == "":
                continue
            else:
                spl = line.split(":")

                if len(spl) != 4:
                    raise Exception(
                        f"Sequence line {i}: expected 4 colon-separated values, found {len(spl)}"
                    )

                res, resolution, inv, rep = spl

                if resolution == "AA":
                    resolution_lookup = "atomistic"
                elif resolution == "UA":
                    resolution_lookup = "united_atom"
                elif resolution == "CG":
                    resolution_lookup = "coarse_grained"
                else:
                    raise Exception(
                        f"Unsupported resolution {resolution}; only AA, UA and CG are currently supported"
                    )

                index.append(
                    monomers.index({"name": res, "resolution": resolution_lookup})
                )

                try:
                    inverted.append(bool(int(inv)))
                except ValueError:
                    raise Exception(
                        f"Sequence line {i}: cannot convert {inv} (Inverted field) to a boolean value"
                    )

                try:
                    rep_int = int(rep)
                except ValueError:
                    raise Exception(
                        f"Sequence line {i}: cannot convert {rep} (Repeats field) to an integer value"
                    )

                if rep_int < 0:
                    raise Exception(
                        f"Sequence line {i}: the Repeats field should contain a positive integer"
                    )

                reps.append(rep_int)

    return Sequence(index, inverted, reps)


def calc_minimal_distance(coords1: np.array, coords2: np.array) -> float:
    """Calculates the minimal distance between all points in a given vector of shape (N,3)
       or the minimal distance between two vectors of shape (N,3) and (M,3). Optimized for speed with
       scipy.spatial.distance.cdist, so that even for large vectors the calculation is fast.

    Args:
        - `coords1 (np.array)`: vector of shape (N,3)
        - `coords2 (np.array)`: vector of shape (M,3)

    Returns:
        - `(float)`: minimal distance between all points in the given vectors
    """

    distance_matrix = cdist(coords1, coords2)

    # set distance of points to themselves to inf
    distance_matrix[distance_matrix < 2 * eps] = np.inf

    return np.min(distance_matrix)


def get_rand_shift(
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
       - `x_min (float)`: Lower bound for the x direction.
       - `x_max (float)`: Upper bound for the x direction.
       - `y_min (float)`: Lower bound for the y direction.
       - `y_max (float)`: Upper bound for the y direction.
       - `z_min (float)`: Lower bound for the z direction.
       - `z_max (float)`: Upper bound for the z direction.
       - `z_bias (float)`: Bias for the z direction.

    Returns:
        - `(np.ndarray)`: 3d vector with shape (3,)
    """
    x = np.random.uniform(x_min, x_max)
    y = np.random.uniform(y_min, y_max)
    z = z_bias * np.random.uniform(z_min, z_max)
    return np.array([x, y, z])


def get_semi_random_walk_shift(
    polypeptide_coordinates: np.array,
    monomer: Monomer,
    trr: float,
    cshift: List[float],
    shift_conf: List[float],
    max_iter: int = 1e3,
) -> np.ndarray:
    """Self-Avoiding Random Walk to prevent infinite forces upon energy minimization.
        Combines RandShift() and MinimalDistance(). Returns k, which can be used to add onto the new monomer.
        Treshshold of trr=1 should be sufficient to prevent infinite forces. If the function fails to find a valid shift, it will return None.

    Args:
        - `polypetide_coordinates (np.ndarray, (3,N))`: Array of all coordinates of the atoms in the previous monomers.
        - `monomer (Monomer)`: Monomer from raccoon.src.data
        - `trr (float)`: Threshold for minimal distance.
        - `cshift (List(float))`: Cartesian shift.
        - `shift_conf (List(float))`: Shift configuration for RandShift()
        - `max_iter (int)`: Maximum number of iterations.

    Returns:
        - `k (np.ndarray)`: 3d vector with shape (3,) if successful, else False.
    """

    monomer_coordinates = monomer.coordinates_to_numpy() + cshift
    minimal_distance = 0
    cnt = 0

    while minimal_distance < trr and cnt < max_iter:
        k = get_rand_shift(*shift_conf)
        updated_monomer_coordinates = monomer_coordinates + k

        minimal_distance = calc_minimal_distance(
            coords1=polypeptide_coordinates, coords2=updated_monomer_coordinates
        )

        cnt += 1

    return k if cnt < max_iter else False


def generate_file(
    monomers: Monomers,
    sequence: Sequence,
    explicit_bonds: bool,
    outpath: str,
    trr: float = 1,
    shift_conf: List[float] = [-1, 1, -1, 1, -1, 1, 1],
    damping_factor: float = 0.5,
    suppress_messages: bool = True,
    cycle_cnt: int = 0,
    max_cycles: int = 100,
):
    """Central function of the modul: adds monomers to a polymer peptide chain and writes it directly to a PDB file.
    If more than 10 cycles are needed to generate a valid structure, the function will raise an error.

    Args:
       - `monomers (Monomers)`: Monomers object.
       - `sequence (Sequence)`: Sequence object.
       - `explicit_bonds (bool)`: If True, explicit bonds are generated.
       - `outpath (str)`: Path to the output file.
       - `trr (float)`: Threshold for minimal distance.
       - `shift_cartesian (list(float))`: Cartesian shift of dimensions x,y,z.
       - `damping_factor (float)`: Damping factor for the shift.
       - `suppress_messages (bool)`: If True, messages are suppressed.
       - `cycle_cnt (int)`: Current cycle count.
       - `max_cycles (int)`: Maximum number of cycles.
    """

    if cycle_cnt >= max_cycles:
        raise Exception(
            f"""Exceeded maximum number of cycles ({max_cycles}) to generate a valid structure.\n
                Please reparameterize the semi random walk. You can increase the treshold `trr`,\n
                reduce the damping factor `damping_factor`, play with the parameters in to generate\n
                the random shift vector `shift_conf`, or increase `max_cycles`."""
        )

    atom_count = 0
    res_count = 0

    with open(outpath, "w") as f:
        cshift = np.zeros(3)
        atom_count = 0

        links = list()
        links_explicit = list()
        coordinates = np.zeros((1, 3))

        if not suppress_messages:
            console = Console()
            console.print("Generating Coordinates")

        for index, inverted, reps in zip(
            sequence.index, sequence.inverted, sequence.reps
        ):
            monomer = monomers[index]

            if inverted:
                monomer = monomer.invert()

            for rep in range(reps):
                shift_conf[6] = float(monomer.atom_count) * damping_factor

                m = get_semi_random_walk_shift(
                    coordinates,
                    monomer,
                    trr=trr,
                    cshift=cshift,
                    shift_conf=shift_conf,
                )

                if not np.all(m):
                    generate_file(
                        monomers,
                        sequence,
                        explicit_bonds,
                        outpath,
                        trr,
                        shift_conf,
                        damping_factor,
                        suppress_messages,
                        cycle_cnt + 1,
                        max_cycles,
                    )

                    # end the current function call and continue with the next one
                    return

                cshift = cshift + m

                updated_monomer = monomer.update(atom_count, cshift)

                coordinates = np.vstack(
                    (coordinates, updated_monomer.coordinates_to_numpy())
                )

                pairs = []

                # Probably not the most efficient method to get all the explicit links, it works however
                # updated_monomer.get_explicit_links() gets a list of all links of every atom in a monomer
                # iterate through that list and create pairs

                if explicit_bonds == True:
                    if not suppress_messages:
                        console.print("Generating Explicit Bonds")

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
                    # formatting coordinates to 3 decimal places
                    x = f"{atom.x:.3f}"
                    y = f"{atom.y:.3f}"
                    z = f"{atom.z:.3f}"

                    f.write(
                        "{:>0}{:<7}{:<5}{:<5}{:<4}{:<3}{:<6}{:<8}{:<8}{:<10}{:<7}{:<14}{}\n".format(
                            "",
                            "ATOM",
                            atom.index,
                            atom.ff_identifier,
                            updated_monomer.name,
                            "A",
                            res_count,
                            x,
                            y,
                            z,
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
        if not suppress_messages:
            console.print("Writing to File")
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

    if not suppress_messages:
        console.print(
            f"Created PDB file with {res_count} residue and {atom_count} atoms to {outpath}.",
            style="bold green",
        )


def close_PDB(fpath: str, atom_count: int):
    """
    Closes a PDB file by adding the END statement and the number of atoms.

    Args:
        - `file_path (str)`: Path to the file.
        - `atom_count (int)`: Number of atoms.

    """

    with open(fpath, "a") as file:
        file.write(
            "{:>0}{:<11}{:<5}{:<5}{:<4}{:<5}{:<5}{:<5}{:<5}{:<5}{:<5}{:<5}{:<5}{}\n".format(
                "", "MASTER", 0, 0, 0, 0, 0, 0, 0, 0, atom_count, 0, atom_count, 0
            )
        )
        file.write("END")
