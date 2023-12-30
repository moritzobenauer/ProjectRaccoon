from unittest import TestCase

from raccoon.src.typing import List, Dict

from raccoon.src.util import MONOMERFILE

from raccoon.src.data import Monomer, Monomers, Atom

import importlib

import json


class TestMonomerJsonClass(TestCase):
    def setUp(self) -> None:
        with open(
            importlib.resources.files("raccoon.src.data") / MONOMERFILE, "r"
        ) as f:
            self.monomers_json = json.load(f)
        return super().setUp()

    def test_resolution(self) -> None:
        """Test if all resolutions are valid."""

        for _, monomer in self.monomers_json.items():
            self.assertIn(
                monomer["resolution"], ["atomistic", "united_atom", "coarse_grained"]
            )

    # test if all resolutions are in atomisitc, united atomes or cg?
    def test_names(self) -> None:
        """Test if all Names are valid."""

        for _, monomer in self.monomers_json.items():
            self.assertEqual(monomer["name"], monomer["name"].upper())
