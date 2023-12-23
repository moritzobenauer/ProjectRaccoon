from typing import List, Tuple, Union, Dict, Optional, Slice, Set
import copy
import ast

import numpy as np


class Monomer:
    """
    Contains the monomer classes for the raccoon model.
    """

    name: str
    """Name of the monomer."""
    resolution: str
    """Resolution of the monomer."""
    atom_count: int
    """Number of atoms in the monomer."""
    atoms: List[Tuple[str, str]]
    """List of atoms in the monomer."""
    link: List[int]  # TODO: rename to links
    """List of atoms that are linked."""
    polymer: bool = False
    """"""
    inverted: bool = False
    """"""

    def __init__(
        self,
        name: str,
        resolution: str,
        atom_count: int,
        atoms: List[List[Union[str, int]]],
        link: List[int],
        polymer: bool = False,
        inverted: bool = False,
    ):
        self.name = name
        self.resolution = resolution
        self.atom_count = atom_count
        self.atoms = atoms
        self.link = link
        self.polymer = polymer
        self.inverted = inverted

    @classmethod
    def from_dict(cls, data: Dict):
        """
        Creates a monomer from a dictionary.

        """

        # convert data to correct types
        data["polymer"] = bool(int(data["polymer"]))
        data["atom_count"] = int(data["atoms"])
        data["link"] = ast.literal_eval(data["link"])
        for x in range(1, data["atoms"] + 1):
            data[x] = ast.literal_eval(data.pop(f"{x}"))
        data["inverted"] = False

        # alle restlichen einträge, deren keys integer sind, in die atoms liste speichern
        data["atoms_list"] = list()
        for key in data.keys():
            if key not in [
                "res",
                "resolution",
                "atoms",
                "link",
                "polymer",
                "inverted",
                "atoms_list",
            ] and isinstance(eval(key), int):
                data["atoms_list"].append([ast.literal_eval(data[key])])

        monomer = cls(
            name=data["res"],
            resolution=data["resolution"],
            atom_count=data["atoms"],
            atoms=data["atoms_list"],
            link=data["link"],
            polymer=data["polymer"],
            inverted=data["inverted"],
        )

        return monomer.invert_amino_acid()

    def invert_amino_acid(self):
        """
        Inverts an amino acid by reversing the link list and changing the inverted flag.

        Returns:
            Monomer: Inverted monomer
        """

        inv_monomer = copy.deepcopy(self)

        inv_monomer.link = inv_monomer.link[::-1]
        inv_monomer.inverted = not inv_monomer.inverted

        return inv_monomer

    def update(self, shift: int, shift_cartesian: List[float]) -> "Monomer":
        """
        Updates the monomer by shifting the atom positions and indicies.

        Args:
            shift (int): Shift value
            shift_cartesian (function): Function to shift cartesian coordinates

        Returns:
            Monomer: Updated monomer
        """

        updated_monomer = copy.deepcopy(self)

        updated_monomer.link = [x + shift for x in updated_monomer.link]

        for atom in updated_monomer.atoms:
            # shifting the atom positions
            for i in range(2, 5):
                atom[i] = np.round(atom[i] + shift_cartesian[i - 3], 2)

            # shifting the atom index
            atom[6] += shift

            # shifting the neighboring atoms
            for i in range(len(atom[5])):
                atom[5][i] += shift

        return updated_monomer

    def __repr__(self):
        return f"Monomer({self.name}, {self.atoms}, {self.link}, {self.polymer})"

    def __str__(self):
        return f"Monomer({self.name}, {self.atoms}, {self.link}, {self.polymer})"

    def __eq__(self, other):
        if isinstance(other, Monomer):
            return self.name == other.name
        else:
            return False

    def __hash__(self):
        return hash(self.name)

    def __getitem__(self, key):
        if isinstance(key, str):
            return self.atoms[self.atoms.index(key)]
        else:
            return self.atoms[key]

    def __len__(self):
        return len(self.atoms)

    def __iter__(self):
        return iter(self.atoms)

    def __contains__(self, item):
        return item in self.atoms


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
        Creates a monomer from a file.
        """

        monomers = []
        with open(fpath) as inpt:
            monomer = dict()
            for line in inpt:
                if line.startswith("#") or line.strip() == "":
                    if monomer:
                        monomers.append(Monomer.from_dict(monomer))
                        monomer = dict()
                    continue
                else:
                    key, value = line.split("=")
                    key = key.strip()
                    value = value.strip()

                    monomer[key] = value

        return cls(monomers)

    def __repr__(self):
        return f"{len(self.monomers)} Monomers"

    def __getitem__(
        self, index: Union[int, List[int], slice]
    ) -> Union[Monomer, List[Monomer]]:
        if isinstance(index, int):
            return self.monomers[index]
        elif isinstance(index, list):
            assert all(isinstance(i, int) for i in index)
            return [self.monomers[i] for i in index]
        elif isinstance(index, slice):
            return self.monomers[index]
        else:
            raise TypeError("Index must be int, list or slice.")

    def __len__(self):
        return len(self.monomers)

    def __iter__(self):
        return iter(self.monomers)

    def __contains__(self, item):
        return item in self.monomers

    def __eq__(self, other: Union[Monomer, Dict]):
        if isinstance(other, Monomer):
            for attribute in self.__dict__:
                if getattr(self, attribute) != getattr(other, attribute):
                    return False
            return True
        elif isinstance(other, Dict):
            for attribute in other.keys():
                if not hasattr(self, attribute):
                    print(f"Monomer does not have attribute {attribute}")
                if getattr(self, attribute) != other[attribute]:
                    return False
            return True

    def __index__(self, other: Union[Monomer, Dict]):
        for monomer in self.monomers:
            if monomer == other:
                return self.monomers.index(monomer)
