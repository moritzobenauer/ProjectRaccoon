import unittest
from pathlib import Path

if __name__ == "__main__":

    root_dir = Path(__file__).resolve().parents[2]

    loader = unittest.TestLoader()
    testSuite = loader.discover(start_dir=root_dir, pattern="test_*.py")
    runner = unittest.TextTestRunner()
    runner.run(testSuite)
