from unittest import TestCase

from raccoon.src.typing import List, Dict

from raccoon.src.util import MONOMERFILE

from raccoon.src.data import Monomer, Monomers
import json

from pathlib import Path


class TestMonomerClass(TestCase):

    """Test the Monomer class."""

    root = Path(__file__).parents[3]

    def setUp(
        self, fpath: str = "raccoon/tests/unit/data/test_monomer_from_dict.json"
    ) -> None:
        """Create a monomer object."""

        with open(Path.joinpath(self.root, fpath), "r") as f:
            self.monomer = Monomer.from_dict(json.load(f))
        return super().setUp()

    def test_prepare_dict(
        self, fpath: str = "raccoon/tests/unit/data/test_prepare_dict_for_monomer.json"
    ) -> None:
        """
        Tests the prepare_dict method with the example PEO Monomer.
        """
        with open(Path.joinpath(self.root, fpath), "r") as f:
            clean_dict = Monomer.prepare_dict(json.load(f))

        self.assertEqual(self.monomer.name, "PEO")

        self.assertEqual(self.monomer.resolution, "united_atom")

        self.assertEqual(self.monomer.atom_count, 3)

        self.assertEqual(self.monomer.link, [3, 1])

        self.assertTrue(self.monomer.polymer)

        self.assertEqual(self.monomer, Monomer.from_dict(clean_dict))

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
            "link": [3, 1],
        }

        false_other = {
            "name": "PEO",
            "resolution": "united_atom",
            "link": [2, 1],
        }

        self.assertEqual(self.monomer.name, other["name"])
        self.assertEqual(self.monomer.resolution, other["resolution"])

        self.assertEqual(self.monomer.link, other["link"])
        self.assertEqual(self.monomer.atom_count, other["atom_count"])

        self.assertEqual(self.monomer, other)
        self.assertNotEqual(self.monomer, false_other)

    def test_invert(
        self, other: str = "raccoon/tests/unit/data/test_invert_monomer.json"
    ) -> None:
        inv_monomer = self.monomer.invert()

        with open(Path.joinpath(self.root, other), "r") as f:
            inv_test_monomer = Monomer.from_dict(json.load(f))

        self.assertEqual(inv_monomer, inv_test_monomer)

        self.assertTrue(inv_monomer.inverted)

    def test_get_atoms_from_bs_file(self) -> None:
        """Tests the get_atoms_from_bs_file method."""
        pass

    def test_create_monomer(self) -> None:
        """Tests the create_monomer method."""
        pass

    def test_update(self) -> None:
        pass


class TestMonomersClass(TestCase):

    """Test the Monomers Class."""

    root = Path(__file__).parents[3]

    def setUp(self) -> None:
        self.monomers = Monomers.from_file()
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
        self.assertEqual(self.monomers.index(self.monomers[1].__dict__), 1)

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
        self.assertTrue(
            peo.name + "_" + peo.resolution in self.monomers.to_dict().keys()
        )

        self.assertTrue(peo in self.monomers)
