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

    def test_update(self) -> None:
        pass

    def test_add_to_file(self, fpath: str = "out"):
        pass


class TestMonomersClass(TestCase):

    """Test the Monomers Class."""

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
