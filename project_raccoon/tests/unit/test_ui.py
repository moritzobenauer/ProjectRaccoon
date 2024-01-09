from unittest import TestCase

from unittest.mock import patch

from project_raccoon.src.data import Sequence
from project_raccoon.src.ui import start_racoon
from project_raccoon.src.ui.user_interface import choose_option, manage_monomers
from project_raccoon.src.functions import get_elements_and_coords_from_pdb


import tempfile
import shutil
import importlib
from pathlib import Path


class TestUserInterface(TestCase):
    def setUp(self) -> None:
        self.seq_file = (
            importlib.resources.files("project_raccoon.tests.unit.data")
            / "seq_FHFHFXG_PEO_GXFHFHF.txt"
        )
        self.seq_file_name = "seq_FHFHFXG_PEO_GXFHFHF.txt"
        self.pdb_file = (
            importlib.resources.files("project_raccoon.tests.unit.data") / "out.pdb"
        )
        self.pdb_file_name = "out.pdb"
        self.xyz_file = (
            importlib.resources.files("project_raccoon.tests.unit.data") / "out.xyz"
        )
        self.xyz_file_name = "out.xyz"

        return super().setUp()

    @patch("project_raccoon.src.ui.user_interface.choose_option", side_effect=["Exit"])
    @patch("builtins.print")
    def test_exit_option(self, *args):
        """Test that the exit option works."""
        result = start_racoon(
            "sequence_file", "out_file", None, True, False, suppress_messages=True
        )

        self.assertIsNone(result)

    @patch(
        "project_raccoon.src.ui.user_interface.choose_option",
        side_effect=["Manage Monomers", "Exit"],
    )
    @patch(
        "project_raccoon.src.ui.user_interface.manage_monomers", side_effect=["Return"]
    )
    @patch("builtins.print")
    def test_return_manage_monomers(self, *args):
        """Test that the manage monomers option works."""
        result = start_racoon(
            "sequence_file", "out_file", None, True, False, suppress_messages=True
        )

        self.assertIsNone(result)

    @patch(
        "project_raccoon.src.ui.user_interface.choose_option",
        side_effect=["Create PDB File", "Exit"],
    )
    @patch("builtins.print")
    def test_create_pdb_file(self, *args):
        """Tests the "Create PDB File" option."""
        with tempfile.TemporaryDirectory() as tmpdir:
            shutil.copy(self.seq_file, tmpdir)

            stmp = Path(tmpdir) / self.seq_file_name

            otmp = Path(tmpdir) / self.pdb_file_name
            result = start_racoon(
                stmp,
                otmp.__str__(),
                None,
                True,
                False,
                suppress_messages=True,
            )

            self.assertIsNone(result)

            self.assertTrue(otmp.exists())

            elements, _ = get_elements_and_coords_from_pdb(otmp)

            self.assertTrue(len(elements) == 400)

    @patch(
        "project_raccoon.src.ui.user_interface.choose_option",
        side_effect=["Check PDB File", "Exit"],
    )
    @patch("builtins.print")
    def test_check_pdb_file(self, *args):
        """Test that the check pdb file option works. Should return None"""
        with tempfile.TemporaryDirectory() as tmpdir:
            shutil.copy(self.pdb_file, tmpdir)

            otmp = Path(tmpdir) / self.pdb_file_name

            result = start_racoon(
                self.seq_file,
                otmp.__str__(),
                None,
                True,
                False,
                suppress_messages=True,
            )

            self.assertIsNone(result)

    @patch(
        "project_raccoon.src.ui.user_interface.choose_option",
        side_effect=["Convert PDB to XYZ File", "Exit"],
    )
    @patch("builtins.print")
    def test_convert_pdb_to_xyz_file(self, *args):
        """Test that the convert pdb to xyz file option works."""

        with tempfile.TemporaryDirectory() as tmpdir:
            shutil.copy(self.pdb_file, tmpdir)

            otmp = Path(tmpdir) / "out.pdb"
            xyz = Path(tmpdir) / "out.xyz"

            result = start_racoon(
                self.seq_file,
                otmp.__str__(),
                None,
                True,
                False,
                suppress_messages=True,
            )

            self.assertIsNone(result)

            self.assertTrue(xyz.exists())

            with open(xyz, "r") as f:
                lines = f.readlines()

            self.assertTrue(len(lines) == 402)
            self.assertEqual(lines[0], "400\n")

        pass

    @patch(
        "project_raccoon.src.ui.user_interface.choose_option",
        side_effect=["Check Minimal Distance", "Exit"],
    )
    @patch("builtins.print")
    def test_check_minimal_distance(self, *args):
        """Test that the check minimal distance option works."""

        with tempfile.TemporaryDirectory() as tmpdir:
            shutil.copy(self.pdb_file, tmpdir)

            otmp = Path(tmpdir) / self.pdb_file_name

            result = start_racoon(
                self.seq_file,
                otmp.__str__(),
                None,
                True,
                False,
                suppress_messages=True,
            )
        self.assertIsNone(result)

    @patch(
        "project_raccoon.src.ui.user_interface.choose_option",
        side_effect=["Manage Monomers", "Exit"],
    )
    @patch(
        "project_raccoon.src.ui.user_interface.manage_monomers",
        side_effect=["Add Monomer", "Return"],
    )
    def test_add_monomer(self, mock_input, mock_input2):
        """Test that the add monomer option works."""
        pass

    @patch(
        "project_raccoon.src.ui.user_interface.choose_option",
        side_effect=["Manage Monomers", "Exit"],
    )
    @patch(
        "project_raccoon.src.ui.user_interface.manage_monomers",
        side_effect=["Delete Monomer", "Return"],
    )
    def test_delete_monomer(self, mock_input, mock_input2):
        """Test that the delete monomer option works."""
        pass

    @patch(
        "project_raccoon.src.ui.user_interface.choose_option",
        side_effect=["Manage Monomers", "Exit"],
    )
    @patch(
        "project_raccoon.src.ui.user_interface.manage_monomers",
        side_effect=["Print Monomers", "Return"],
    )
    def test_print_monomers(self, mock_input, mock_input2):
        """Test that the print monomers option works."""
        pass

    @patch(
        "project_raccoon.src.ui.user_interface.choose_option",
        side_effect=["Manage Monomers", "Exit"],
    )
    @patch(
        "project_raccoon.src.ui.user_interface.manage_monomers",
        side_effect=["Export JSON Monomer File", "Return"],
    )
    def test_export_json_monomer_file(self, mock_input, mock_input2):
        """Test that the export json monomer file option works."""
        pass
