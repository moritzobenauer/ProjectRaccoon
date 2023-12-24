from unittest import TestCase

from raccoon.src.typing import List, Dict

from raccoon.src.util import MONOMERFILE

from raccoon.src.data import Monomer, Monomers
import json


class TestMonomerClass(TestCase):

    """Test the Monomer class."""

    def setUp(
        self, fpath: str = "raccoon/tests/unit/test_monomer_from_dict.json"
    ) -> None:
        """Create a monomer object."""
        with open(fpath, "r") as f:
            dict = json.load(input)
            self.monomer = Monomer.from_dict(dict)
        return super().setUp()

    def test_prepare_dict(
        self, fpath: str = "raccoon/tests/unit/test_prepare_dict_for_monomer.json"
    ) -> None:
        """
        Tests the prepare_dict method with the example PEO Monomer.
        """
        with open(fpath, "r") as f:
            clean_dict = json.load(input)

        self.assertEqual(self.monomer.name, "PEO")

        self.assertEqual(self.monomer.resolution, "united_atom")

        self.assertEqual(self.monomer.atom_count, 3)

        self.assertEqual(self.monomer.link, [3, 1])

        self.assertTrue(self.monomer.polymer)

        self.assertEqual(self.monomer, Monomer.from_dict(clean_dict))

        pass

    def test_monomer(self) -> None:
        """Tests the containing datatypes."""
        self.assertIsInstance(self.monomer, Monomer)

        self.assertIsInstance(self.monomer.name, str)

        self.assertIsInstance(self.monomer.resolution, str)

        self.assertIsInstance(self.monomer.link, List[int])

        self.assertIsInstance(self.monomer.atom_count, int)

        self.assertIsInstance(
            self.monomer.atoms, List[str, str, float, float, float, List[int], int]
        )

        self.assertIsInstance(self.monomer.polymer, bool)

        self.assertIsInstance(self.monomer.inverted, bool)

        self.assertFalse(self.monomer.inverted)

    def test_equal(self) -> None:
        self.assertTrue(self.monomer == self.monomer)

        other = {
            "name": "PEO",
            "resolution": "united_atom",
            "link": [3, 1],
        }

        false_other = other = {
            "name": "PEO",
            "resolution": "united_atom",
            "link": [2, 1],
        }

        self.assertTrue(self.monomer == other)
        self.assertFalse(self.monomer == false_other)

    def test_invert(
        self, other: str = "raccoon/tests/unit/test_invert_monomer.json"
    ) -> None:
        inv_monomer = self.monomer.invert()

        with open(self.input_file, "r") as f:
            monomer = json.load(input)
            inv_test_monomer = Monomer.from_dict(monomer)

        assert inv_monomer == inv_test_monomer

        self.assertTrue(inv_monomer.inverted)

    def test_update(self, other: str) -> None:
        pass

    def test_add_to_file(self, fpath: str = "out"):
        pass


class TestMonomersClass(TestCase):

    """Test the Monomers Class."""

    def setUp(self) -> None:
        self.monomers = Monomers.from_file()
        return super().setUp()

    def test_get_item(self) -> None:
        pass

    def test_indexing(self) -> None:
        pass
