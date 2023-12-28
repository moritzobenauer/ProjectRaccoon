from biopandas.pdb import PandasPdb

ppdb = PandasPdb().fetch_pdb("3eiy")
import py3Dmol
import pandas as pd
from math import sqrt
from itertools import combinations
from tqdm import tqdm
import os

from sys import float_info

from pathlib import Path


def clear_terminal():
    os.system("cls" if os.name == "nt" else "clear")


def PDBtoXYZ(fpath: str):
    root = Path(__file__).parents[3]

    ppdb.read_pdb(os.path.join(root, fpath))
    ppdb.df["ATOM"]
    columns_to_include = ["element_symbol", "x_coord", "y_coord", "z_coord"]
    n_atoms = len(ppdb.df["ATOM"]) + 1

    fpath = fpath.split(".")[0] + ".xyz"
    with open(os.path.join(root, fpath), "w") as file:
        file.write(str(n_atoms) + "\n")
        file.write("XYZ file generated with Project RACCOON 2023" + "\n")
        for index, row in ppdb.df["ATOM"][columns_to_include].iterrows():
            line = " ".join(map(str, row))
            if index == (n_atoms - 2):
                file.write(line)
            else:
                file.write(line + "\n")
    print(f"PDB file was converted to xyz file {fpath}.")


def CheckPDB(input):
    ppdb.read_pdb(f"{input}.pdb")
    print(ppdb.df["ATOM"])


# This function is only available in a notebook format


def Visualize(input):
    with open(f"{input}.xyz") as ifile:
        system = "".join([x for x in ifile])
    view = py3Dmol.view(width=400, height=300)
    view.addModelsAsFrames(system)
    view.setStyle({"model": -1}, {"sphere": {"color": "spectrum"}})
    view.zoomTo()
    view.show()


def Distances(input):
    ppdb.read_pdb(input)
    ppdb.df["ATOM"]
    columns_to_include = ["element_symbol", "x_coord", "y_coord", "z_coord"]
    for index, row in ppdb.df["ATOM"][columns_to_include].iterrows():
        print(row[1], row[2], row[3])
        print(index + 1)


def calculate_distance(point1, point2):
    return sqrt(
        (point1[0] - point2[0]) ** 2
        + (point1[1] - point2[1]) ** 2
        + (point1[2] - point2[2]) ** 2
    )


def CheckMinimalDistance(fpath: str):
    """
    Checks the minimal distance between two atoms in a pdb file.
    """
    root = Path(__file__).parents[3]
    ppdb.read_pdb(os.path.join(root, fpath))
    columns_to_include = ["element_symbol", "x_coord", "y_coord", "z_coord"]

    # initialize with max float value
    min_distances = float_info.max

    progress_bar = tqdm(
        total=0.5 * len(ppdb.df["ATOM"]) ** 2, desc="Calculating minimal distance"
    )

    for (index1, row1), (index2, row2) in combinations(
        ppdb.df["ATOM"][columns_to_include].iterrows(), 2
    ):
        distance = sqrt(
            (row1.iloc[1] - row2.iloc[1]) ** 2
            + (row1.iloc[2] - row2.iloc[2]) ** 2
            + (row1.iloc[3] - row2.iloc[3]) ** 2
        )

        if min_distances > distance:
            min_distances = distance

        progress_bar.update(1)

    progress_bar.close()
    print(f"Minimal distance: {min_distances:.4f} Å")
