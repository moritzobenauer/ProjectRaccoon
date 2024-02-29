from py3Dmol import view as py3Dmol
import os
import numpy as np

from ..typing import Tuple, List
from .standard import calc_minimal_distance


def get_elements_and_coords_from_pdb(fpath: str) -> Tuple[List[str], np.ndarray]:
    """Returns the elements and coordinates from a pdb file.

    Args:
       - `fpath (str)`: Path to pdb file.

    Returns:
        - `(Tuple[List[str], np.ndarray])`: Elements and coordinates.
    """
    coords = list()
    elements = list()

    with open(fpath, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                atom = line.split(" ")[1:]
                atom = [l.strip() for l in atom if l.strip() != ""]
                coords.append([float(x) for x in atom[5:8]])
                elements.append(atom[-1])

    return elements, np.array(coords)


def get_links_from_pdb(fpath: str) -> np.ndarray:
    """Returns the links from a pdb file.

    Args:
        - `fpath (str)`: Path to pdb file.

    Returns:
        - `(np.ndarray)`: Links.
    """
    links = list()

    with open(fpath, "r") as f:
        for line in f:
            if line.startswith("CONECT"):
                link = line.split(" ")[1:]
                link = [l for l in link if l.strip() != ""]
                links.append([int(l) for l in link])

    return np.array(links)


def pdb_to_xyz(fpath: str, suppress_messages=False) -> None:
    elements, coords = get_elements_and_coords_from_pdb(fpath)

    out_path = fpath.split(".")[0] + ".xyz"

    with open(out_path, "w") as f:
        f.write(str(len(elements)) + "\n")
        f.write("XYZ file generated with Project RACCOON 2023" + "\n")

        for element, coord in zip(elements, coords):
            f.write(f"{element}\t{coord[0]:.2f}\t{coord[1]:.2f}\t{coord[2]:.2f}\n")

    if not suppress_messages:
        print(f"PDB file was converted to xyz file {out_path}.")


def check_pdb_file(fpath: str, suppress_messages: bool = True) -> None:
    """Checks a pdb file for errors with the biopandas package.
       If the file is valid, it prints the atom dataframe.

    Args:
        fpath (str): Path to pdb file.
    """

    from biopandas.pdb import PandasPdb

    ppdb = PandasPdb().fetch_pdb("3eiy")
    ppdb.read_pdb(fpath)

    if not suppress_messages:
        print(ppdb.df["ATOM"])


def visualize_pdb_file(fpath: str) -> None:
    """Visualizes a pdb file with py3Dmol. Only available in a notebook format.

    Args:
        - `fpath (str)`: Path to pdb file.
    """
    view = py3Dmol(width=800, height=500, viewergrid=(1, 1))
    view.addModel(open(f"{fpath}", "r").read(), "pdb")
    view.setStyle({"sphere": {}}, viewer=(0, 0))
    view.zoomTo()
    view.show()


def check_minimal_distance(fpath: str) -> None:
    "Calculates the minimal distance from a given pdb file"
    _, coords = get_elements_and_coords_from_pdb(fpath)

    min_dist = calc_minimal_distance(coords, coords)

    print(f"Minimal distance: {min_dist:.4f} Å")
