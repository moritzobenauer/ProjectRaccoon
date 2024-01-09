from ..typing import List, Dict, Union, Optional
from .structs import Atom

from ..util import MONOMERFILE
import copy
import ast

import importlib

import numpy as np


class Monomer:
    """
    Contains the monomer class for the raccoon model.
    """

    name: str
    """Name of the monomer."""
    resolution: str
    """Resolution of the monomer."""
    atom_count: int
    """Number of atoms in the monomer."""
    atoms: List[Atom]
    """
    List of atoms in the monomer.
    contains the following information:
    [ff_identifier : str , element: str, x: float, y : float, z : float, neighbors : List[int], number of the atom : int]
    """
    link: List[int]  # TODO: ranme to links
    """Index of the linked atoms, C- and N-Terminus."""
    polymer: bool = False
    """"""
    inverted: bool = False
    """"""

    def __init__(
        self,
        name: str,
        resolution: str,
        atom_count: int,
        atoms: List[Atom],
        link: List[int],
        polymer: bool,
        inverted: bool,
    ):
        self.name = name
        self.resolution = resolution
        self.atom_count = atom_count
        self.atoms = atoms
        self.link = link
        self.polymer = polymer
        self.inverted = inverted

    def coordinates_to_numpy(self):
        """Returns a numpy array of the cartesian coordinates of all atoms."""
        return np.array([[atom.x, atom.y, atom.z] for atom in self.atoms])

    def get_explicit_links(self):
        """Returns a list of lists of all explicit links within the monomer."""
        return [atom.neighbours for atom in self.atoms]

    @classmethod
    def prepare_dict(cls, data: Dict):
        """Prepares a dictionary for the creation of a monomer, which is created by
            the `Monomer.from_file()` function, which uses the `monomers.dat` file.

        Args:
            - `data (Dict)`: Dictionary to prepare.

        Returns:
            - `data (Dict)`: Prepared dictionary.
        """
        # convert data to correct types
        data["name"] = data["res"]
        data["polymer"] = bool(int(data["polymer"]))
        data["atom_count"] = int(data["atoms"])
        data["link"] = ast.literal_eval(data["link"])

        # alle restlichen einträge, deren keys integer sind, in die atoms liste speichern
        data["atoms"] = list()
        for key in data.keys():
            if key not in [
                "res",
                "name",
                "resolution",
                "atoms",
                "atom_count",
                "link",
                "polymer",
            ] and isinstance(eval(key), int):
                atom_data = ast.literal_eval(data[key])
                data["atoms"].append(atom_data)

        if not "inverted" in data.keys():
            data["inverted"] = False

        return data

    @classmethod
    def from_dict(cls, data: Dict):
        """
        Creates a monomer from a dictionary.

        Args:
            - `data (Dict)`: Cleaned dictionary to create the monomer from.

        Returns:
            Monomer: Monomer object.

        """

        monomer = cls(
            name=data["name"],
            resolution=data["resolution"],
            atom_count=data["atom_count"],
            atoms=[Atom(*atom) for atom in data["atoms"]],
            link=data["link"],
            polymer=data["polymer"],
            inverted=data["inverted"],
        )

        return monomer

    def invert(self) -> "Monomer":
        """
        Inverts an amino acid by reversing the link list and changing the inverted flag.

        Returns:
            - `inv_monomer (Monomer)`: Inverted monomer
        """
        if self.polymer:
            raise ValueError(
                f"Inversion not possible for {self.name} (polymer building block)."
            )

        inv_monomer = copy.deepcopy(self)

        inv_monomer.link = inv_monomer.link[::-1]
        inv_monomer.inverted = not inv_monomer.inverted

        return inv_monomer

    def update(self, shift: int, shift_cartesian: List[float]) -> "Monomer":
        """
        Updates the monomer by shifting the atom positions and indicies.

        Args:
            - `shift (int)`: Shift value
            - `shift_cartesian (function)`: Function to shift cartesian coordinates

        Returns:
            - `updated_monomer (Monomer)`: Updated monomer
        """

        updated_monomer = copy.deepcopy(self)

        updated_monomer.link = [x + shift for x in updated_monomer.link]

        for atom in updated_monomer.atoms:
            # shifting the atom positions
            atom.x = atom.x + shift_cartesian[0]
            atom.y = atom.y + shift_cartesian[1]
            atom.z = atom.z + shift_cartesian[2]

            # shifting the atom index
            atom.index += shift

            # shifting the neighboring atoms
            for idx, _ in enumerate(atom.neighbours):
                atom.neighbours[idx] += shift

        return updated_monomer

    @classmethod
    def create_monomer(
        cls,
        name: str,
        resolution: str,
        polymer: bool,
        link: List[int],
        atoms: List,
        ff_identifiers: List[str],
    ) -> "Monomer":
        """
        Creates a monomer.
        Args:
            - `name (str)`: Name of the monomer.
            - `resolution (str)`: Resolution of the monomer.
            - `polymer (bool)`: Polymer flag of the monomer.
            - `link (List[int])`: List of atoms that are linked.
            - `atoms (List)`: List of atoms in the monomer from the function Monomer.get_atoms_from_bs_file.
            - `ff_identifier (List[str])`: List of atom names that are used for the force field.

        Returns:
            - `monomer (Monomer)`: Monomer object.
        """

        for atom, ff_identifier in zip(atoms, ff_identifiers):
            atom.insert(0, ff_identifier)

        atoms = [Atom(*atom) for atom in atoms]

        atom_count = len(atoms)

        inverted = False

        return cls(
            name=name,
            resolution=resolution,
            atom_count=atom_count,
            atoms=atoms,
            link=link,
            polymer=polymer,
            inverted=inverted,
        )

    @staticmethod
    def get_atoms_from_bs_file(fpath: str) -> List:
        """Reads a bs file and returns a list of atoms, that are necessary to create a monomer.
           To complete the information the force field identifiers are needed and are added in
           the Monomers.create_monomer function.

        Args:
            - `fpath (str)`: File path to the bs file.

        Returns:
            - `atoms (List)`: List of atoms, which contains the following information: [element: str, x: float, y : float, z : float, neighbors : List[int], number of the atom : int]

        """

        with open(fpath, "r") as f:
            lines = f.readlines()

        atoms = list()
        for idx, line in enumerate(lines[2:], start=1):
            line = line.split()
            atoms.append(
                [
                    line[0],
                    float(line[1]),
                    float(line[2]),
                    float(line[3]),
                    [int(n) for n in line[4:]],
                    idx,
                ]
            )

        return atoms

    def __repr__(self):
        return f"Monomer({self.name}, resolution {self.resolution}, # atoms {self.atom_count}, Polymer:{self.polymer}, inv {self.inverted})"

    def __eq__(self, other: Union["Monomer", Dict]):
        """
        Defines the equality of two monomers by comparing all attributes or an monomer and a dictionary, that contains some of the attributes.
        With this method, the '==' operator is implemented. For the monomers it is important to know, that the monomers are only compared by their name and resolution.

        Args:
            - `other (Union[Monomer, Dict])`: Monomer or dictionary to compare with.

        Returns:
            - `bool`: True if equal, False if not equal.
        """
        if isinstance(other, Monomer):
            for attribute in self.to_dict().keys():
                if attribute in [
                    "name",
                    "resolution",
                    "atom_count",
                    "inverted",
                    "polymer",
                ] and getattr(self, attribute) != getattr(other, attribute):
                    return False
            return True
        elif isinstance(other, Dict):
            for attribute in other.keys():
                if not hasattr(self, attribute):
                    print(f"monomer has no attribute {attribute}")
                    return False
                if attribute in ["atoms", "link"]:
                    continue
                if getattr(self, attribute) != other[attribute]:
                    return False
            return True

    def to_dict(self):
        """Returns a dictionary of the monomer."""
        return {
            "name": self.name,
            "resolution": self.resolution,
            "atom_count": self.atom_count,
            "atoms": [atom.to_list() for atom in self.atoms],
            "link": self.link,
            "polymer": self.polymer,
            "inverted": self.inverted,
        }

    def __hash__(self):
        return hash(self.name)

    def __len__(self):
        return len(self.atoms)

    def __iter__(self):
        return iter(self.atoms)


class Monomers:
    """
    Contains the monomer classes for the raccoon model.
    """

    monomers: List[Monomer]
    """List of monomers."""

    def __init__(self, monomers: List[Monomer]):
        self.monomers = monomers

    @classmethod
    def from_file(cls, fpath: str):
        """
        Creates a list of monomers from a `.dat` file. Is not longer supported and will be removed in the future.

        Args:
            - `fpath (str, optional)`: Path to the monomer file. Defaults to MONOMERFILE.

        Returns:
            - ` monomers (Monomers)`: Monomers object.
        """

        monomers = []

        with open(fpath, "r") as f:
            monomer = dict()

            lines = f.readlines()
            fsize = len(lines)
            for lnr, line in enumerate(lines):
                if line.startswith("#") or not line.strip():
                    if monomer:
                        monomer = Monomer.prepare_dict(monomer)
                        monomers.append(Monomer.from_dict(monomer))
                        monomer = dict()
                    continue
                else:
                    key, value = line.split("=")
                    monomer[key.strip()] = value.strip()

                    if lnr == fsize - 1:
                        # TODO: gefällt mir nicht, andere lösung finden
                        monomer = Monomer.prepare_dict(monomer)
                        monomers.append(Monomer.from_dict(monomer))

        return cls(monomers)

    @classmethod
    def from_json(cls, fpath: Optional[str] = None) -> "Monomers":
        """
        Creates a list of monomers from a json file. If no path is given the default monomers file is used,
        which is located in the project_raccoon package under project_raccoon/src/data/monomers.json.

        Args:
            - `fpath (Optional[str])`: Path to the monomers file. Defaults to None.

        Returns:
            - `monomers (Monomers)`: Monomers object.
        """
        import json

        if fpath is not None:
            with open(fpath, "r") as f:
                data = json.load(f)
        else:
            files = importlib.resources.files("project_raccoon.src.data")
            with open(files / MONOMERFILE, "r") as f:
                data = json.load(f)

        monomers = []

        for _, monomer in data.items():
            monomers.append(Monomer.from_dict(monomer))

        return cls(monomers)

    def add_monomer(self, monomer: Monomer, save: Optional[bool] = False) -> None:
        """
        Adds a monomer to the monomers list, if it is not yet in the list. This is done,
        by checking if another monomer with the same name and resolution is already in the list.
        After that, it is necessary to update the monomers file by calling the Monomers.to_json function.

        Args:
            - `monomer (Monomer)`: Monomer to add.
            - `save (Optional[bool], optional)`: Save the monomer to the monomers file. Defaults to False.
        """

        if {"name": monomer.name, "resolution": monomer.resolution} in self.monomers:
            print("Monomer already in list")
            return

        self.monomers.append(monomer)

        if save:
            self.to_json(MONOMERFILE)

    def remove_monomer(self, monomer: Monomer, save: Optional[bool] = False) -> None:
        """Removes a monomer from the monomers list. After that, it is necessary to update the monomers file by calling the Monomers.to_json function.

        Args:
            - `monomer (Monomer)`: Monomer to remove.
            - `save (Optional[bool], optional)`: Save the monomer to the monomers file. Defaults to False.
        """

        self.monomers.remove(monomer)

        if save:
            self.to_json()

    def __repr__(self):
        return f"{len(self.monomers)} Monomers"

    def __getitem__(
        self, index: Union[int, List[int], slice]
    ) -> Union[Monomer, "Monomers"]:
        if isinstance(index, int):
            return self.monomers[index]
        elif isinstance(index, list):
            assert all(isinstance(i, int) for i in index)
            return Monomers([self.monomers[i] for i in index])
        elif isinstance(index, slice):
            return Monomers(self.monomers[index])
        else:
            raise TypeError("Index must be int, list or slice.")

    def to_dict(self):
        """Returns a dictionary of the monomers."""
        return {
            monomer.name + "_" + monomer.resolution: monomer.to_dict()
            for monomer in self.monomers
        }

    def to_json(self, fpath: Optional[str] = None, indent: Optional[int] = 2):
        """Saves the monomers to a json file. If no path is given the default monomers file is used,
           which is located in the project_raccoon package under project_raccoon/src/data/monomers.json.

        Args:
            - `fpath (Optional[str])`: Path to the monomers file. Defaults to None.
            - `indent (Optional[int], optional)`: Indentation of the json file. Defaults to 2.

        """
        import json

        if fpath is not None:
            with open(fpath, "w") as f:
                json.dump(self.to_dict(), f, indent=indent)
        else:
            files = importlib.resources.files("project_raccoon.src.data")
            with open(files / MONOMERFILE, "w") as f:
                json.dump(self.to_dict(), f, indent=indent)

    def __len__(self):
        return len(self.monomers)

    def __iter__(self):
        return iter(self.monomers)

    def __contains__(self, item):
        return item in self.monomers

    def index(self, monomer: Monomer) -> int:
        """
        Return the index of an monomer in the monomers.
        """
        if monomer in self.monomers:
            return self.monomers.index(monomer)
        else:
            raise ValueError(
                f"{monomer.name} in {monomer.resolution} monomer.resolution not in list"
            )

    def __sizeof__(self) -> int:
        return len(self.monomers)

    def __len__(self):
        return len(self.monomers)
