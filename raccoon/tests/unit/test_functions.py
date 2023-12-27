from unittest import TestCase

from raccoon.src.functions import generate_file, generate_sequence

from raccoon.src.typing import List, Dict, Tuple, NamedTuple

from raccoon.src.data import Monomer, Monomers, Sequence, Atom
from collections import namedtuple

from pathlib import Path

import numpy as np


class TestFunctions(TestCase):
    def setUp(self) -> None:
        self.root = Path(__file__).parents[3]
        self.monomers = Monomers.from_file()
        return super().setUp()

    def test_generate_sequence(
        self,
        spath: str = "raccoon/tests/unit/data/seq_FHFHFXG_PEO_GXFHFHF.txt",
    ) -> None:
        s1 = generate_sequence(self.monomers, spath)

        index = [6, 3, 4, 3, 4, 3, 9, 7, 1, 7, 9, 3, 4, 3, 4, 3, 6]
        inverted = [True] * len(index)
        reps = [1, 1, 1, 1, 1, 1, 1, 1, 50, 1, 1, 1, 1, 1, 1, 1, 1]

        s2 = Sequence(index, inverted, reps)

        self.assertIsInstance(s1, Sequence)

        self.assertIsInstance(s1.index, List)
        self.assertIsInstance(s1.index[0], int)

        self.assertIsInstance(s1.inverted, List)
        self.assertIsInstance(s1.inverted[0], bool)

        self.assertIsInstance(s1.reps, List)
        self.assertIsInstance(s1.reps[0], int)
        self.assertEqual(s1, s2)

    def test_generate_file(self):
        pass
