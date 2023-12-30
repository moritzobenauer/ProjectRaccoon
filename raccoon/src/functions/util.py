import py3Dmol
import pandas as pd
from math import sqrt
from itertools import combinations
from tqdm import tqdm
import os

from scipy.spatial.distance import cdist
import numpy as np


from sys import float_info

from pathlib import Path


def clear_terminal():
    os.system("cls" if os.name == "nt" else "clear")


def get_elements_and_coords_from_pdb(fpath: str):
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


def pdb_to_xyz(fpath: str) -> None:
    elements, coords = get_elements_and_coords_from_pdb(fpath)

    out_path = fpath.split(".")[0] + ".xyz"

    with open(out_path, "w") as f:
        f.write(str(len(elements)) + "\n")
        f.write("XYZ file generated with Project RACCOON 2023" + "\n")

        for element, coord in zip(elements, coords):
            f.write(f"{element}\t{coord[0]:.2f}\t{coord[1]:.2f}\t{coord[2]:.2f}\n")

        print(f"PDB file was converted to xyz file {out_path}.")


def CheckPDB(input):
    from biopandas.pdb import PandasPdb

    ppdb = PandasPdb().fetch_pdb("3eiy")
    ppdb.read_pdb(f"{input}.pdb")
    print(ppdb.df["ATOM"])


# This function is only available in a notebook format


def Visualize(input):
    view = py3Dmol.view(width=800, height=500, viewergrid=(1, 1))
    view.addModel(open(f"{input}", "r").read(), "pdb")
    view.setStyle({"sphere": {}}, viewer=(0, 0))
    view.zoomTo()
    view.show()


def calculate_distance(point1, point2):
    return sqrt(
        (point1[0] - point2[0]) ** 2
        + (point1[1] - point2[1]) ** 2
        + (point1[2] - point2[2]) ** 2
    )


def calc_minimal_distance(coords: np.array) -> float:
    """calculates the minimal distance between all points in a given vector"""

    distance_matrix = cdist(coords, coords)

    np.fill_diagonal(distance_matrix, np.inf)

    return np.min(distance_matrix)


def CheckMinimalDistance(fpath: str):
    "calculates the minimal distance from a given pdb file"
    _, coords = get_elements_and_coords_from_pdb(fpath)

    coords = np.array(coords)

    min_dist = calc_minimal_distance(coords)

    print(f"Minimal distance: {min_dist:.4f} Å")
