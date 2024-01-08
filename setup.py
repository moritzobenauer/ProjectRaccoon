from setuptools import setup, find_packages
from pathlib import Path

long_description = Path(Path(__file__).parent / "README.md").read_text()
requirements = Path(Path(__file__).parent / "requirements.txt").read_text().splitlines()

setup(
    name="raccoon",
    version="1.0.0",
    description="Project RACCOON is a Python tool designed for the generation of PDB (Protein Data Bank) files for polymer peptide conjugates, polypeptides, and polymers in a building block fashion.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Moritz L. Obenauer, Kai N. Spauszus.",
    url="https://github.com/moritzobenauer/ProjectRaccoon",
    packages=find_packages(),
    package_data={
        "raccoon": ["src/data/*.json", "requirements.txt", "README.md", "LICENSE"]
    },
    test_suite="raccoon/tests",
    entry_points={"console_scripts": ["raccoon=raccoon.__main__:main"]},
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
        "Natural Language :: English",
    ],
)
