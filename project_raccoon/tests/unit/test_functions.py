from contextlib import AbstractContextManager
from typing import Any
from unittest import TestCase

from project_raccoon.src.functions import (
    generate_file,
    generate_sequence,
    calc_minimal_distance,
    get_elements_and_coords_from_pdb,
    get_links_from_pdb,
    pdb_to_xyz,
)
from project_raccoon.src.functions.standard import get_semi_random_walk_shift

from project_raccoon.src.typing import List, Dict, Tuple, NamedTuple

from project_raccoon.src.data import Monomer, Monomers, Sequence, Atom
from collections import namedtuple

from pathlib import Path

import numpy as np

import tempfile
import shutil
import importlib


class TestFunctions(TestCase):
    def setUp(self) -> None:
        self.seq_file_name = "seq_FHFHFXG_PEO_GXFHFHF.txt"
        self.seq_file = (
            importlib.resources.files("project_raccoon.tests.unit.data")
            / self.seq_file_name
        )
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
    ) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            shutil.copy(self.seq_file, Path(tmpdir) / self.seq_file_name)
            seq = generate_sequence(self.monomers, Path(tmpdir) / self.seq_file_name)

        self.assertIsInstance(seq, Sequence)

        self.assertIsInstance(seq.index, List)
        self.assertIsInstance(seq.index[0], int)

        self.assertIsInstance(seq.inverted, List)
        self.assertIsInstance(seq.inverted[0], bool)

        self.assertIsInstance(seq.reps, List)
        self.assertIsInstance(seq.reps[0], int)

        self.assertEqual(len(seq.index), len(self.seq.index))

        self.assertEqual(seq.inverted, self.seq.inverted)
        self.assertEqual(seq.reps, self.seq.reps)

    def test_srw_rs(self) -> None:
        """Test the semi random walk shift with a random shift vector"""

        trr = 1
        shift_conf = [-1, 1, -1, 1, -1, 1, 1]
        damping_factor = 0.5

        atom_count = 0
        res_count = 0

        # cartesian shifts
        cshift = np.zeros(3)
        atom_count = 0

        coordinates = np.zeros((1, 3))

        sequence = self.seq
        for index, _, reps in zip(sequence.index, sequence.inverted, sequence.reps):
            monomer = self.monomers[index]

            for rep in range(reps):
                shift_conf[6] = float(monomer.atom_count) * damping_factor
                m = get_semi_random_walk_shift(
                    coordinates,
                    monomer,
                    trr=trr,
                    cshift=cshift,
                    shift_conf=shift_conf,
                )

                cshift += m

                updated_monomer = monomer.update(atom_count, cshift)

                coordinates = np.vstack(
                    (coordinates, updated_monomer.coordinates_to_numpy())
                )

                atom_count += updated_monomer.atom_count
                res_count += 1

        min_dist = calc_minimal_distance(coordinates, coordinates)

        self.assertTrue(min_dist >= trr)

    def test_calc_minimal_distance(self) -> None:
        """Tests the minimal distance function"""

        coordinates = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 3], [4, 5, 6]])

        # test N,3 with N,3
        min_dist = calc_minimal_distance(coordinates, coordinates)
        self.assertAlmostEqual(min_dist, np.sqrt(3))

        # test N,3 with 1,3
        min_dist = calc_minimal_distance(coordinates[:-1], coordinates[-1:])
        self.assertAlmostEqual(min_dist, 4.69041575982343)

    def test_srw_nrs(self) -> None:
        """Test the semi random walk shift with a non random shift vector"""
        # 2. test SRW
        #   - random shift minmal distance > threshold
        #   - test non random shift cumsum np.arange(len(atoms))
        pass

    def test_explicit_bonds(self) -> None:
        """Tests the explicit bond function"""
        pass

    def test_get_elements_and_coords_from_pdb(self) -> None:
        pass

    def test_get_links_from_pdb(self) -> None:
        pass
