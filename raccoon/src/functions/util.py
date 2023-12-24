from biopandas.pdb import PandasPdb

ppdb = PandasPdb().fetch_pdb("3eiy")
import py3Dmol
import pandas as pd
import numpy as np
from itertools import combinations
from tqdm import tqdm


def PDBtoXYZ(input):
    ppdb.read_pdb(f"{input}.pdb")
    ppdb.df["ATOM"]
    columns_to_include = ["element_symbol", "x_coord", "y_coord", "z_coord"]
    n_atoms = len(ppdb.df["ATOM"]) + 1
    with open(f"{input}.xyz", "w") as file:
        file.write(str(n_atoms) + "\n")
        file.write("XYZ file generated with Project RACCOON 2023" + "\n")
        for index, row in ppdb.df["ATOM"][columns_to_include].iterrows():
            line = " ".join(map(str, row))
            if index == (n_atoms - 2):
                file.write(line)
            else:
                file.write(line + "\n")
    print("PDB file was converted to xyz file.")


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
    return np.sqrt(
        (point1[0] - point2[0]) ** 2
        + (point1[1] - point2[1]) ** 2
        + (point1[2] - point2[2]) ** 2
    )


def CheckMinimalDistance(input):
    ppdb.read_pdb(f"{input}.pdb")
    columns_to_include = ["element_symbol", "x_coord", "y_coord", "z_coord"]
    distances = []

    total_iter = (ppdb.df["ATOM"][columns_to_include].shape[0] ** 2) / 2
    progress_bar = tqdm(total=total_iter, desc="Processing")

    for (index1, row1), (index2, row2) in combinations(
        ppdb.df["ATOM"][columns_to_include].iterrows(), 2
    ):
        point1 = (float(row1[1]), float(row1[2]), float(row1[3]))
        point2 = (float(row2[1]), float(row2[2]), float(row2[3]))
        distance = calculate_distance(point1, point2)
        distances.append(distance)
        progress_bar.update(1)

    progress_bar.close()
    print(min(distances))
