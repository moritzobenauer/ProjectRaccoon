![raccoon](/screenshots/asset1.png)

# Project RACCOON

 [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black) ![JGU Mainz](https://img.shields.io/badge/JGU%20Mainz%20-%20red.svg) [![Coverage](https://codecov.io/gh/moritzobenauer/ProjectRaccoon/coverage.svg?branch=master)](https://codecov.io/gh/moritzobenauer/ProjectRaccoon?branch=master)


**Automated construction of atomistic and coarse-grained models in the PDB format for linear polymer peptide conjugates.**

## General Purpose & Scope
 Project RACCOON (*Rapid Automated Construction of Conjugates using the Optimized OPLS Input*) is a Python tool designed for the generation of PDB (Protein Data Bank) files for polymer peptide conjugates, polypeptides, and polymers in a building block fashion. It allows for the easy addition of new monomers, incorporation of polypeptide and polymer sequences in text form, and outputs PDB and XYZ files.
The tool was specifically developed for the OPLS force field but is adaptable to any other force field. The tool's primary scope is to provide a straightforward method for creating starting structures for molecular dynamics simulations. It caters to the creation of linear polymer peptide conjugates with facing polypeptide strands. 
The tool emphasizes user-friendliness, making it accessible to chemists and physicists, not just limited to theoretical chemists with extensive background knowledge. Additionally, users can modify the code for their specific projects. Early developments started in January 2023. The software was developed at the KOMET and Besenius' research groups at the Johannes Gutenberg University Mainz.

### Supported building blocks

>[!NOTE]
>We encourage users to share their updated *monomers.json* file. Please submit a pull request, open an issue, or contact the developers.

| Residue Name | AA | UA | CG |
| :---         | :- | :- | :- |
| PEO          | :white_check_mark: | :white_check_mark: | :soon: |
| GLY          | :white_check_mark: | :soon: | :soon: |
| PHE          | :white_check_mark: | :white_check_mark: | :soon: |
| HIS          | :white_check_mark: | :white_check_mark: | :soon: |
| ACE          | :white_check_mark: | :white_check_mark: | :x: |
| NME          | :soon: | :soon: | :x: |
| DUM          |  |  | :white_check_mark: |
| CSX          | :white_check_mark: | :white_check_mark: | :x: |
| LNK          | :white_check_mark: | :white_check_mark: | :x: |

For non-standard amino acids **CSX** and **LNK**, please refer to the following paper by [pending, 20xx](https://www.ak-besenius.chemie.uni-mainz.de/). Coarse-grained models have not been implemented yet.



## How to use Project RACCOON

### Installing for Standard Usage

The easiest way to use the Project RACCOON software is to install it via pip

```bash
pip install project_raccoon
```

If you want to have access to the gromax test functions, you can alternatively clone the github repository


```bash
git clone https://github.com/moritzobenauer/ProjectRaccoon.git
```

and install the modul via

```bash
pip install -e .
```

To start the command-line interface run

```bash
project_raccoon
```

### Importing new Monomers

Open a graphical software to create three-dimensional molecular models and obtain the desired monomer structure. It is important to note that only the atoms or beads of the repeating units are used and not the entire molecule. The structure must be saved in .bs format, as the position of the atoms or beads and the neighboring atoms or beads are specified here. Run the command 
```
python import_monomer -f newmonomer.bs -o monomers.dat
```
to add the new monomer unit to your library. Please enter the name of the new residue and some properties and give every atom a unique identifier according to the force field you plan to use. For amino acid residues, it is essential to specify the C- and N-terminus of the building block.

![Bild1](/screenshots/import1.png)
![Bild2](/screenshots/import2.png)



### Creating PDB Files from Sequence
The basis for creating a PDB file is a seq.txt file containing the block sequence. The monomers are specified line by line with their name and resolution (this may be necessary if the same block is present in the library file with a different resolution). In addition, whether it is an inverted building block (this occurs with C2 symmetrical polymer peptide conjugates) and the number of repeating units must be specified. In this way, classical polymers with a fixed number of repeating units of a monomer can also be represented quickly.
```
HIS:AA:0:10
PHE:AA:0:10
HIS:AA:0:10
```
The example shown here corresponds to a polypeptide in which ten histidine units are followed by ten phenylalanines and then ten histidines again. All monomers are shown in atomistic resolution and are not inverted.

| Residue Name | Resolution   | Inverted | Repeats |
| :---         | :---         | :---     | :---    |
| String       | {AA, UA, CG} | {0,1}    | Integer |

The program is started with the following command.
```
python raccoon.py -s seq.txt -o resultingfile -m monomers.dat
```
To create the PDB file, select the first item that appears in the selection. The PDB file is then generated automatically according to the selected parameters. For visualization or other purposes, creating an XYZ file from the PDB file may be desired. This XYZ file no longer contains any binding information. For this purpose, the option *Convert PDB to XYZ file* can be selected in the menu. A corresponding XYZ file with the same file name as the PDB file is created.


![Bild3](/screenshots/raccoon_main.png)


| Parameter | Default       |   Description |
| :---      | :---          | :---          |
| -s        | seq.txt       | Filename with extension for the sequence.                            |
| -o        | out           | Filename without extension for the resulting pdb file.               |
| -m        | monomers.dat  | Filename with extension for the monomer library.                     |
| -e        | False         | Write all explicit bonds. It can be useful for nonstandard residues. |
| -r        | True          | Remove duplicate bonds. Bond a --> b is equivalent to b --> a.       |

### Check PDB Files

The successful creation of the PDB file can then be checked. The *Check PDB* file option is selected for this purpose. If the PDB file appears in tabular form in the terminal, it is correct and will be read correctly by all standard programs. Individual atoms or beads may be arranged too close to each other. In the worst case, this can lead to problems with energy minimization or simulations. To check that no two atoms or beads are too close to each other, the *Check Minimal Distance* option can be selected. The output contains the smallest distance between two atoms or beads.
## Examples

The software was developed specifically for the telechelic polymer peptide conjugates, according to Otter et. al. 2018. The peptide sequence FHFHFXG-PEO-GXFHFHF (with X: 6-aminohexanoic acid) is an example. A corresponding seq.txt file can be found in the examples folder. A graphical representation of the PDB file is shown in the image below. However, it should be noted that the arrangement of the monomers here is in a parabolic form (with no loss of generality). 

![Bild4](/screenshots/raccoon_export.png)

After a brief energy minimization (performed here with GROMACS and the OPLS force field), the physically meaningful structure shown in the following image is obtained. The PDB file created with Project RACCOON is an excellent starting point for such molecular dynamics simulations. It is important to recognize here that different building blocks with atomistic and united atom resolution have been combined with each other. 

![Bild5](/screenshots/raccoon_em.png)


## Limitations

> [!WARNING]
> * The generated output is only a starting structure for molecular dynamics simulations.
> * Users must ensure atom/bead naming consistency for matching force field data.

