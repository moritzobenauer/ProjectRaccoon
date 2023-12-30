from unittest import TestCase

# 1. teste allgemeine eigenschaften
#   - anzahl der atome, gleich denen der sequenz sum atom_count*rep for atom in seq
#   - file appropriately closed
#   - sorted
#   - links


class TestPdbFile(TestCase):

    """
    läuft wenn in main check pdf_file aufgerufen wird
    """

    def setUp(self) -> None:
        return super().setUp()

    def test_is_sorted(self):
        pass

    def test_closing(self):
        pass

    def test_validate_pdb_file(self):
        pass
