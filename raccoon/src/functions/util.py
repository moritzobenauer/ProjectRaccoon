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


## PDB2XYZ can cause problems
## Use PDB format for visualization


# def PDBtoXYZ(fpath: str):
#    from biopandas.pdb import PandasPdb
#
#    ppdb = PandasPdb().fetch_pdb("3eiy")
#    root = Path(__file__).parents[3]
#
#    ppdb.read_pdb(os.path.join(root, fpath))
#    ppdb.df["ATOM"]
#    columns_to_include = ["element_symbol", "x_coord", "y_coord", "z_coord"]
#    n_atoms = len(ppdb.df["ATOM"]) + 1
#
#    fpath = fpath.split(".")[0] + ".xyz"
#    with open(os.path.join(root, fpath), "w") as file:
#        file.write(str(n_atoms) + "\n")
#        file.write("XYZ file generated with Project RACCOON 2023" + "\n")
#        for index, row in ppdb.df["ATOM"][columns_to_include].iterrows():
#            line = " ".join(map(str, row))
#            if index == (n_atoms - 2):
#                file.write(line)
#            else:
#                file.write(line + "\n")
#    print(f"PDB file was converted to xyz file {fpath}.")


def pdb_to_xyz(fpath: str) -> None:
    coords = list()
    elements = list()

    with open(fpath, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                atom = line.split(" ")[1:]
                atom = [l.strip() for l in atom if l.strip() != ""]
                coords.append([float(x) for x in atom[5:8]])
                elements.append(atom[-1])
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


def Distances(input):
    from biopandas.pdb import PandasPdb

    ppdb = PandasPdb().fetch_pdb("3eiy")
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

    from biopandas.pdb import PandasPdb

    ppdb = PandasPdb().fetch_pdb("3eiy")
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
