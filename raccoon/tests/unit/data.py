from unittest import TestCase

from raccoon.src.typing import List, float, int, slice, Dict

from raccoon.src.util import MONOMERFILE

from raccoon.src.data import Monomer, Monomers
import json


class TestMonomerClass(TestCase):

    """Test the Monomer class."""

    input_file: str = "raccoon/tests/test_data/test_monomer.json"

    def setUp(self) -> None:
        """Create a monomer object."""
        with open(self.input_file, "r") as f:
            dict = json.load(input)
            self.monomer = Monomer.from_dict(dict)
        return super().setUp()

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
        self.assertIsNot(self.monomer.inverted)

    def test_equal(self, other: str) -> None:
        pass

    def test_invert(self, other: str) -> None:
        inv_monomer = self.monomer.invert()
        with open(self.input_file, "r") as f:
            monomer = json.load(input)

            inv_test_monomer = Monomer.from_dict(monomer)

            assert inv_monomer == inv_test_monomer

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
