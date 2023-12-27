from ..typing import List, Dict, Union, Optional

from ..util import MONOMERFILE
import copy
import ast

import questionary

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
    atoms: List
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
        atoms: List,
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

    @classmethod
    def prepare_dict(cls, data: Dict):
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
                data["atoms"].append(ast.literal_eval(data[key]))

        if not "inverted" in data.keys():
            data["inverted"] = False

        return data

    @classmethod
    def from_dict(cls, data: Dict):
        """
        Creates a monomer from a dictionary.

        """

        monomer = cls(
            name=data["name"],
            resolution=data["resolution"],
            atom_count=data["atom_count"],
            atoms=data["atoms"],
            link=data["link"],
            polymer=data["polymer"],
            inverted=data["inverted"],
        )

        return monomer

    def invert(self) -> "Monomer":
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

    def add_to_file(self, fpath: str = MONOMERFILE):
        """
        Adds a monomer to a given file, default the monomers.dat file.

        Args:
            fpath (str, optional): Path to the monomer file. Defaults to MONOMERFILE.

        """

        monomer = [
            "\n",
            "res=" + self.name + "\n",
            "resolution=" + self.resolution + "\n",
            "polymer=" + str(int(self.polymer)) + "\n",
            "atoms=" + str(self.atom_count) + "\n",
            "links=" + str(self.link) + "\n",
        ]

        with open(fpath, "a") as f:
            f.writelines(monomer)

            for idx, atom in enumerate(self.atoms):
                line = f"""{idx}=["{atom[0]}", "{atom[1]}", {atom[3]}, {atom[4]}, {atom[5]}, {atom[2]}, {idx}] \n"""
                f.write(line)

        print(f"Monomer added to {fpath}")

    @classmethod
    def create_monomer(
        cls,
        fpath: str,
        monomer_name: Optional[str] = None,
        monomer_resolution: Optional[str] = None,
        monomer_polymer: Optional[bool] = None,
        monomer_link: Optional[List[int]] = None,
        ff_identifiers: Optional[List[int]] = None,
    ) -> "Monomer":
        """
        Creates a monomer. The function can be automated by passing the optional arguments and/or can be used interactively.
        Args:
            fpath (str): Path to the monomer file.
            monomer_name (Optional[str], optional): Name of the monomer. Defaults to None.
            monomer_resolution (Optional[str], optional): Resolution of the monomer. Defaults to None.
            monomer_polymer (Optional[bool], optional): Polymer flag of the monomer. Defaults to None.
            monomer_link (Optional[List[int]], optional): List of atoms that are linked. Defaults to None.
            ff_identifier (Optional[List[int]], optional): List of atom indices that are used for the force field. Defaults to None.
        """

        monomer = Monomer()

        if not monomer_name:
            monomer_name = questionary.text("Enter the name of the monomer").ask()

        if not monomer_resolution:
            monomer_resolution = questionary.select(
                "Choose resolution",
                choices=["atomistic", "united_atom", "coarse_grained"],
            ).ask()

        if not monomer_polymer:
            monomer_polymer = questionary.confirm("Is this a polymer?").ask()

        atoms = Monomer.get_atoms_from_bs_file(fpath)

        if not ff_identifiers:
            for atom in atoms:
                ff_identifier = questionary.text(
                    f"Enter a force field Identifier for {atom[0]}' with the Number {atom[-1]}: "
                ).ask()
                atom.insert(0, ff_identifier)

        atom_count = len(atoms)

        if not monomer_link:
            linkC = questionary.text(f"Choose C-Terminus (1-{len(atoms)})").ask()
            linkN = questionary.text(f"Choose N-Terminus (1-{len(atoms)})").ask()
            monomer_link = [int(linkC), int(linkN)]

        link = monomer_link

        inverted = False

        return cls(
            name=monomer_name,
            resolution=monomer_resolution,
            atom_count=atom_count,
            atoms=atoms,
            link=link,
            polymer=monomer_polymer,
            inverted=inverted,
        )

    @staticmethod
    def get_atoms_from_bs_file(fpath: str) -> List:
        """Reads a bs file and returns a list of atoms, that are necessary to create a monomer.
           To complete the information the force field identifiers are needed and are added in
           the Monomers.create_monomer function.

        Args:
            fpath (str): File path to the bs file.

        Returns:
            List: List of atoms, which contains the following information: [element: str, x: float, y : float, z : float, neighbors : List[int], number of the atom : int]

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
        return f"Monomer({self.name}, resolution {self.resolution}, # Atome {self.atom_count},Polymer:{self.polymer}, inv {self.inverted})"

    def __eq__(self, other: Union["Monomer", Dict]):
        """
        Defines the equality of two monomers by comparing all attributes or an monomer and a dictionary, that contains some of the attributes.
        With this method, the '==' operator is implemented. For the monomers it is important to know, that the monomers are only compared by their name and resolution.

        Args:
            other (Union[Monomer, Dict]): Monomer or dictionary to compare with.

        Returns:
            bool: True if equal, False if not equal.
        """
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

    def to_dict(self):
        return self.__dict__

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
    def from_file(cls, fpath: str = MONOMERFILE):
        """
        Creates a list of monomers from a `.dat` file. Is not longer supported and will be removed in the future.

        Args:
            fpath (str, optional): Path to the monomer file. Defaults to MONOMERFILE.

        Returns:
            Monomers: Monomers object.
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
    def from_json(cls, fpath: str) -> "Monomers":
        """
        Creates a list of monomers from a json file.
        """
        import json

        with open(fpath, "r") as f:
            data = json.load(f)

        monomers = []

        for name, monomer in data.items():
            monomers.append(Monomer.from_dict(monomer))

        return cls(monomers)

    def add_monomer(self, monomer: Monomer, save: Optional[bool] = False) -> None:
        """
        Adds a monomer to the monomers list, if it is not yet in the list. This is done,
        by checking if another monomer with the same name and resolution is already in the list.
        After that, it is necessary to update the monomers file by calling the Monomers.to_json function.

        Args:
            monomer (Monomer): Monomer to add.
            save (Optional[bool], optional): Save the monomer to the monomers file. Defaults to False.
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
            monomer (Monomer): Monomer to remove.
            save (Optional[bool], optional): Save the monomer to the monomers file. Defaults to False.
        """

        self.monomers.remove(monomer)

        if save:
            self.to_json(MONOMERFILE)

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
        return {
            monomer.name + "_" + monomer.resolution: monomer.to_dict()
            for monomer in self.monomers
        }

    def to_json(self, fpath: str, indent: int = 2):
        import json

        with open(fpath, "w") as f:
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
            raise ValueError(f"{monomer} not in list")

    def __sizeof__(self) -> int:
        return len(self.monomers)

    # define len(monomers) to len(self.monomers)
    def __len__(self):
        return len(self.monomers)
