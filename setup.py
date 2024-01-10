from setuptools import setup, find_packages

from pathlib import Path

parent = Path(__file__).parent.resolve()

long_description = (parent / "README.md").read_text()

setup(
    name="project_raccoon",
    version="1.0.1",
    description="Project project_raccoon is a Python tool designed for the generation of PDB (Protein Data Bank) files for polymer peptide conjugates, polypeptides, and polymers in a building block fashion.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Moritz L. Obenauer, Kai N. Spauszus.",
    url="https://github.com/moritzobenauer/ProjectRaccoon",
    packages=find_packages(),
    package_data={
        "project_raccoon": ["src/data/monomers.json"],
    },
    test_suite="project_raccoon/tests",
    entry_points={"console_scripts": ["project_raccoon=project_raccoon.__main__:main"]},
    python_requires='>=3.11',
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
