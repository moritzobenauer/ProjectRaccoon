from unittest import TestCase

from raccoon.src.functions import (
    generate_file,
    generate_sequence,
    calc_minimal_distance,
)
from raccoon.src.functions.standard import SemiRandomWalk

from raccoon.src.typing import List, Dict, Tuple, NamedTuple

from raccoon.src.data import Monomer, Monomers, Sequence, Atom
from collections import namedtuple

from pathlib import Path

import numpy as np


class TestFunctions(TestCase):
    def setUp(self) -> None:
        self.root = Path(__file__).parents[3]
        self.monomers = Monomers.from_json()
        index = [6, 3, 4, 3, 4, 3, 9, 7, 1, 7, 9, 3, 4, 3, 4, 3, 6]
        inverted = [
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            True,
            True,
            True,
            True,
            True,
            True,
            True,
            False,
        ]
        reps = [1, 1, 1, 1, 1, 1, 1, 1, 50, 1, 1, 1, 1, 1, 1, 1, 1]

        self.seq = Sequence(index, inverted, reps)

        return super().setUp()

    def test_generate_sequence(
        self,
        spath: str = "raccoon/tests/unit/data/seq_FHFHFXG_PEO_GXFHFHF.txt",
    ) -> None:
        seq = generate_sequence(self.monomers, spath)

        self.assertIsInstance(seq, Sequence)

        self.assertIsInstance(seq.index, List)
        self.assertIsInstance(seq.index[0], int)

        self.assertIsInstance(seq.inverted, List)
        self.assertIsInstance(seq.inverted[0], bool)

        self.assertIsInstance(seq.reps, List)
        self.assertIsInstance(seq.reps[0], int)
        self.assertEqual(seq, self.seq)

    def test_srw_rs(self) -> None:
        """Test the semi random walk shift with a random shift vector"""

        trr = 1
        shift_cartesian = [-1, 1, -1, 1, -1, 1, 1]
        damping_factor = 0.5

        atom_count = 0
        res_count = 0

        # cartesian shifts
        cshifts = np.zeros(3)
        atom_count = 0

        coordinates = np.zeros((1, 3))

        sequence = self.seq
        for index, _, reps in zip(sequence.index, sequence.inverted, sequence.reps):
            monomer = self.monomers[index]

            for rep in range(reps):
                shift_cartesian[6] = float(monomer.atom_count) * damping_factor
                m = SemiRandomWalk(coordinates, monomer, trr=trr, shift=shift_cartesian)

                cshifts += m

                updated_monomer = monomer.update(atom_count, cshifts)

                new_coordinates = updated_monomer.coordinates_to_numpy()
                coordinates = np.vstack((coordinates, new_coordinates))

                atom_count += updated_monomer.atom_count
                res_count += 1

        min_dist = calc_minimal_distance(coordinates, coordinates)

        self.assertEqual(min_dist, trr)
        self.assertTrue(min_dist >= trr)

    def test_srw_nrs(self) -> None:
        """Test the semi random walk shift with a non random shift vector"""
        pass
        # 2. test SRW
        #   - random shift minmal distance > threshold
        #   - test non random shift cumsum np.arange(len(atoms))

    def test_explicit_bonds(self) -> None:
        """Tests the explicit bond function"""
        pass
