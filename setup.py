from setuptools import setup, find_packages
from pathlib import Path

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
    packages=find_packages(),
    package_data={"raccoon": ["src/data/*.json"]},
    test_suite="raccoon/tests",
    entry_points={"console_scripts": ["raccoon=raccoon.__main__:main"]},
    install_requires=[
        "biopandas==0.4.1",
        "numpy==1.26.3",
        "pandas==2.1.4",
        "py3Dmol==2.0.4",
        "questionary==2.0.1",
        "rich==13.7.0",
        "scipy==1.11.4",
        "setuptools==68.2.2",
        "tqdm==4.66.1",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
        "Natural Language :: English",
    ],
)
