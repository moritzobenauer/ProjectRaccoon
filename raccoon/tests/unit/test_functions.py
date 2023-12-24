from unittest import TestCase

from raccoon.src.functions import generate_file, generate_sequence

from raccoon.src.typing import Monomers, List, Dict, Tuple, NamedTuple

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
        f1path: str = "raccoon/tests/unit/data/seq_FHFHFXG_PEO_GXFHFHF.txt",
    ) -> None:
        with open(Path.joinpath(self.root, f1path), "r") as f:
            s1 = generate_sequence(self.monomers, f1path)

        s2 = namedtuple("sequence", ["index", "inverted", "reps"])
        s2.index = [6, 3, 4, 3, 4, 3, 9, 7, 1, 7, 9, 3, 4, 3, 4, 3, 6]
        s2.inverted = [True] * len(s2.index)
        s2.reps = [1, 1, 1, 1, 1, 1, 1, 1, 50, 1, 1, 1, 1, 1, 1, 1, 1, 1]

        self.assertIsInstance(s1, NamedTuple)

        self.assertIsInstance(s1.index, List)
        self.assertIsInstance(s1.index[0], int)

        self.assertIsInstance(s1.inverted, List)
        self.assertIsInstance(s1.inverted[0], bool)

        self.assertIsInstance(s1.reps, List)
        self.assertIsInstance(s1.reps[0], int)
        self.assertEqual(s2, s2)

    def test_generate_file(self):
        pass
