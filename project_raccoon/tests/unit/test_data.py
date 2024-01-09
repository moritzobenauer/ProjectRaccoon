from unittest import TestCase

from project_raccoon.src.typing import List, Dict

from project_raccoon.src.util import MONOMERFILE

from project_raccoon.src.data import Monomer, Monomers, Atom
import numpy as np

from pathlib import Path


class TestMonomerClass(TestCase):

    """Test the Monomer class."""

    root = Path(__file__).parents[3]

    def assertListAlmostEqual(self, list1: List, list2: List, places=5) -> None:
        """Asserts that two lists are almost equal."""
        self.assertEqual(len(list1), len(list2))
        for a, b in zip(list1, list2):
            self.assertAlmostEqual(a, b, places=places)

    def setUp(self) -> None:
        """Create a monomer object."""
        monomer_dict = {
            "name": "PEO",
            "resolution": "united_atom",
            "atom_count": 3,
            "atoms": [
                ["CA", "C", 1.18, -0.011, -0.29, [2, 0], 1],
                ["CB", "C", 1.09, -0.01, 1.19, [1, 3], 2],
                ["O", "O", -0.08, 0.5, 1.44, [4, 8], 3],
            ],
            "link": [3, 1],
            "polymer": True,
            "inverted": False,
        }

        self.monomer = Monomer.from_dict(monomer_dict)
        return super().setUp()

    def test_prepare_dict(self) -> None:
        """
        Tests the prepare_dict method with the example PEO Monomer.
        """
        clean_dict = Monomer.prepare_dict(
            {
                "res": "PEO",
                "resolution": "united_atom",
                "polymer": "1",
                "atoms": "3",
                "link": "[3, 1]",
                "1": '[ "CA", "C", 1.18, -0.011, -0.29, [ 2, 0 ], 1 ]',
                "2": '[ "CB", "C", 1.09, -0.010,  1.19, [ 1, 3 ], 2 ]',
                "3": '[  "O", "O",-0.08,  0.500,  1.44, [ 4, 8 ], 3 ]',
            }
        )

        monomer = Monomer.from_dict(clean_dict)

        self.assertEqual(monomer.name, "PEO")

        self.assertEqual(monomer.resolution, "united_atom")

        self.assertEqual(monomer.atom_count, 3)

        self.assertEqual(monomer.link, [3, 1])

        self.assertIsInstance(monomer.polymer, bool)
        self.assertTrue(monomer.polymer)

        self.assertEqual(monomer.atom_count, len(monomer.atoms))

        self.assertIsInstance(monomer.atoms, List)
        self.assertIsInstance(monomer.atoms[0], Atom)

        self.assertFalse(monomer.inverted)

    def test_monomer(self) -> None:
        """Tests the containing datatypes."""
        self.assertIsInstance(self.monomer, Monomer)

        self.assertIsInstance(self.monomer.name, str)

        self.assertIsInstance(self.monomer.resolution, str)

        self.assertIsInstance(self.monomer.link, List)
        self.assertIsInstance(self.monomer.link[0], int)

        self.assertIsInstance(self.monomer.atom_count, int)

        self.assertIsInstance(self.monomer.atoms, List)

        self.assertTrue(len(self.monomer.atoms) == self.monomer.atom_count)

        self.assertIsInstance(self.monomer.polymer, bool)

        self.assertIsInstance(self.monomer.inverted, bool)

        self.assertFalse(self.monomer.inverted)

    def test_equal(self) -> None:
        self.assertEqual(self.monomer, self.monomer)

        other = {
            "name": "PEO",
            "resolution": "united_atom",
            "atom_count": 3,
            "inverted": False,
            "polymer": True,
            "link": [3, 1],
        }

        false_other_1 = {
            "name": "PEO",
            "resolution": "united_atom",
            "atom_count": 4,
            "inverted": False,
            "polymer": True,
        }

        false_other_2 = {
            "name": "PEO",
            "resolution": "united_atom",
            "atom_count": 3,
            "inverted": True,
            "polymer": True,
        }

        false_other_3 = {
            "name": "PEO",
            "resolution": "united_atom",
            "atom_count": 3,
            "inverted": False,
            "polymer": False,
        }

        false_other_4 = {
            "name": "test_PEO",
            "resolution": "united_atom",
            "atom_count": 3,
            "inverted": False,
            "polymer": False,
        }

        false_other_5 = {
            "name": "test_PEO",
            "resolution": "united_atom",
            "atom_count": 4,
            "invert": False,
            "polymer": True,
        }

        self.assertEqual(self.monomer.name, other["name"])
        self.assertEqual(self.monomer.resolution, other["resolution"])

        self.assertEqual(self.monomer.link, other["link"])
        self.assertEqual(self.monomer.atom_count, other["atom_count"])

        self.assertEqual(self.monomer, other)
        self.assertNotEqual(self.monomer, false_other_1)  # test resolution
        self.assertNotEqual(self.monomer, false_other_2)  # test invert
        self.assertNotEqual(self.monomer, false_other_3)  # test polymer
        self.assertNotEqual(self.monomer, false_other_4)  # test name
        self.assertNotEqual(self.monomer, false_other_5)  # test atom_count

    def test_invert(self) -> None:
        with self.assertRaises(ValueError) as context:
            self.monomer.invert()  # polymer blocks (monomer.polymer == True) cannot be inverted, peo is polymer

        monomer = Monomer.from_dict(
            {
                "name": "LNK",
                "resolution": "united_atom",
                "atom_count": 10,
                "atoms": [
                    ["NB", "N", -0.777, 0.372, -0.458, [2, 7], 1],
                    ["CB", "C", 0.453, 1.051, -0.071, [1, 3], 2],
                    ["C", "C", 1.706, 0.173, -0.014, [2, 4, 5], 3],
                    ["NA", "N", 2.859, 0.936, -0.15, [3, 6, 10], 4],
                    ["OXT", "O", 1.744, -1.035, 0.2, [3], 5],
                    ["HNA", "H", 2.724, 1.932, -0.056, [4], 6],
                    ["HNB", "H", -0.851, -0.488, 0.101, [1], 7],
                    ["CD", "C", 4.159, 0.34, 0.057, [4, 9], 8],
                    ["CE", "C", 5.297, 1.328, -0.18, [8, 10], 9],
                    ["O", "O", 6.559, 0.707, 0.047, [9], 10],
                ],
                "link": [10, 1],
                "polymer": False,
                "inverted": False,
            }
        )

        inv_monomer = monomer.invert()

        self.assertEqual(inv_monomer.link, monomer.link[::-1])  # links are swtiched

        self.assertTrue(inv_monomer.inverted)  # inverted == True

    def test_get_atoms_from_bs_file(self) -> None:
        """Tests the get_atoms_from_bs_file method."""
        pass

    def test_create_monomer(self) -> None:
        """Tests the create_monomer method."""

        name = "test"
        resolution = "test_resolution"
        polymer = False
        link = [1, 2]
        atoms = [
            ["C", 0, 0, 0, [1, 2], 1],
            ["C", 0, 0, 0, [2, 3], 2],
            ["C", 0, 0, 0, [3, 1], 3],
        ]

        ff_identifier = "C1", "C2", "C3"

        monomer = Monomer.create_monomer(
            name=name,
            resolution=resolution,
            polymer=polymer,
            link=link,
            atoms=atoms,
            ff_identifiers=ff_identifier,
        )

        self.assertIsInstance(monomer, Monomer)

    def test_update(self) -> None:
        shift = 5
        shift_cartesian = [1.0, 1.0, 1.0]

        updated_monomer = self.monomer.update(
            shift=shift, shift_cartesian=shift_cartesian
        )

        self.assertIsInstance(updated_monomer, Monomer)
        self.assertIsInstance(updated_monomer.atoms, List)
        self.assertIsInstance(updated_monomer.atoms[0], Atom)

        # test updating atom index
        self.assertEqual(
            self.monomer.atoms[0].index + shift, updated_monomer.atoms[0].index
        )
        self.assertListAlmostEqual(
            [atom.index + shift for atom in self.monomer.atoms],
            [atom.index for atom in updated_monomer.atoms],
        )

        # test updating cartesian coordinates

        x = [atom.x + shift_cartesian[0] for atom in self.monomer.atoms]
        x_ = [atom.x for atom in updated_monomer.atoms]

        self.assertIsInstance(x, List)
        self.assertIsInstance(x[0], float)
        self.assertIsInstance(x_, List)
        self.assertIsInstance(x_[0], float)

        self.assertAlmostEqual(
            self.monomer.atoms[0].x + shift_cartesian[0], updated_monomer.atoms[0].x
        )

        self.assertListAlmostEqual(x, x_)

        # test updating neighbours
        self.assertListAlmostEqual(
            [
                [neighbour + shift for neighbour in atom.neighbours]
                for atom in self.monomer.atoms
            ],
            [atom.neighbours for atom in updated_monomer],
        )

    def test_coordinates_to_numpy(self) -> np.ndarray:
        """Tests the coordinates_to_numpy method."""
        coordinates = self.monomer.coordinates_to_numpy()
        self.assertIsInstance(coordinates, np.ndarray)
        self.assertEqual(coordinates.shape, (self.monomer.atom_count, 3))

    def test_get_explicit_links(self) -> List[List[int]]:
        """Tests the get_explicit_links method."""
        links = self.monomer.get_explicit_links()

        self.assertEqual(len(links), len(self.monomer.atoms))
        self.assertIsInstance(links, List)
        self.assertIsInstance(links[0], List)
        self.assertIsInstance(links[0][0], int)

    def test_repr(self) -> None:
        """Tests the __repr__ method."""
        self.assertIsInstance(self.monomer.__repr__(), str)

    def test_len(self) -> None:
        """Tests the __len__ method."""
        self.assertEqual(len(self.monomer), self.monomer.atom_count)


class TestMonomersClass(TestCase):

    """Test the Monomers Class."""

    root = Path(__file__).parents[3]

    def setUp(self) -> None:
        self.monomers = Monomers.from_json()
        return super().setUp()

    def test_get_item(self) -> None:
        """Tests the get_item method with index, list of indices and a slice."""
        self.assertIsInstance(self.monomers[0], Monomer)

        monomers_list = self.monomers[0:2]

        self.assertIsInstance(monomers_list, Monomers)
        self.assertIsInstance(monomers_list[0], Monomer)

        monomers_list = self.monomers[[2, 3, 5]]

        self.assertIsInstance(monomers_list, Monomers)
        self.assertIsInstance(monomers_list[0], Monomer)

    def test_indexing(self) -> None:
        """Tests the indexing of the monomers with an monomer and a dict."""
        self.assertEqual(self.monomers.index(self.monomers[0]), 0)
        self.assertEqual(self.monomers.index(self.monomers[2].to_dict()), 2)

    def test_to_dict(self) -> None:
        """
        Tests the to_dict method. Preliminary test with the original `monomers.dat` file,
         which will be removed in the future. Only the `monomers.json` file will be used, tested and updated.
        """

        monomers = Monomers.from_file("monomers.dat")

        self.assertIsInstance(self.monomers.to_dict(), Dict)

        self.assertEqual(len(monomers), len(self.monomers))

        self.assertEqual(self.monomers.to_dict(), monomers.to_dict())

    def test_add_monomer(self) -> None:
        """Tests the add_monomer method."""

        name = "PEO"
        resolution = "test_resolution"
        link = [1, 2]
        atom_count = 3
        atoms = [
            ["C1", "C", 0, 0, 0, [1, 2], 1],
            ["C2", "C", 0, 0, 0, [2, 3], 2],
            ["C3", "C", 0, 0, 0, [3, 1], 3],
        ]
        atoms = [Atom(*atom) for atom in atoms]
        polymer = True
        inverted = False

        peo = Monomer(
            name=name,
            resolution=resolution,
            link=link,
            atom_count=atom_count,
            atoms=atoms,
            polymer=polymer,
            inverted=inverted,
        )

        self.monomers.add_monomer(peo)

        # identifier for monomer is name + resolution
        d = self.monomers.to_dict()
        self.assertTrue(peo.name + "_" + peo.resolution in d.keys())

        self.assertTrue(peo in self.monomers)
