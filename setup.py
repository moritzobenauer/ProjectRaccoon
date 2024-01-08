from setuptools import setup
from pathlib import Path


HERE = Path(__file__).parent

with open(HERE / "requirements.txt") as f:
    required = f.read().splitlines()

setup(
    name="raccoon",
    version="1.0.0",
    description="Project RACCOON\
                 (Rapid Automated Construction of Conjugates using the Optimized OPLS Input)\
                 is a Python tool designed for the generation of PDB (Protein Data Bank) files\
                 for polymer peptide conjugates, polypeptides, and polymers in a building block fashion.\
                 It allows for the easy addition of new monomers, incorporation of polypeptide and\
                 polymer sequences in text form, and outputs PDB and XYZ files.",
    author="Moritz L. Obenauer, Kai N. Spauszus.",
    url="https://github.com/moritzobenauer/ProjectRaccoon",
    packages=["raccoon"],
    package_data={"raccoon": ["src/data/*.json"]},
    test_suite="raccoon/tests",
    entry_points={"console_scripts": ["raccoon=raccoon.__main__:main"]},
    install_requires=required,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
)
