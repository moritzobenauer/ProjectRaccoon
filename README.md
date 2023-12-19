<div>
<img style="display: flex" src="/screenshots/raccoon_logo_round.png" width="100" height="100">
<h1 style="float: right">Project RACCOON</h1> 
</div>
**Automated construction of atomistic and coarse-grained models in the PDB format for linear polymer peptide conjugates.**

## General Purpose & Scope
 <div align="justify">Project RACCOON (<i>Rapid Automated Construction of Conjugates using the Optimized OPLS Input</i>) is a Python tool designed for the generation of PDB (Protein Data Bank) files for polymer peptide conjugates, polypeptides, and polymers in a building block fashion. It allows for the easy addition of new monomers, incorporation of polypeptide and polymer sequences in text form, and outputs PDB and XYZ files.
The tool was specifically developed for the OPLS force field but is adaptable to any other force field. The tool's primary scope is to provide a straightforward method for creating starting structures for molecular dynamics simulations. It caters to the creation of linear polymer peptide conjugates with facing polypeptide strands. 
The tool emphasizes user-friendliness, making it accessible to chemists and physicists, not just limited to theoretical chemists with extensive background knowledge. Additionally, users can modify the code for their specific projects. Early developments started in January 2023. The software was developed at the KOMET and Besenius' research groups at the Johannes Gutenberg University Mainz.</div>

## How to use RACCOON

### Importing new Monomers

Open a graphical software to create three-dimensional molecular models and obtain the desired monomer structure. It is important to note that only the atoms or beads of the repeating units are used and not the entire molecule. The structure must be saved in .bs format, as the position of the atoms or beads and the neighboring atoms or beads are specified here. Run the command 
```
python import_monomer -f newmonomer.bs -o monomers.dat
```
to add the new monomer unit to your library. Please enter the name of the new residue and some properties and give every atom a unique identifier according to the force field you plan to use. For amino acid residues, it is essential to specify the C- and N-terminus of the building block.

<div>
<img style="display: inline-block" src="/screenshots/import1.png" width="600" height="400">
<img style="display: inline-block" src="/screenshots/import2.png" width="600" height="400">
</div>


### Creating PDB Files from Sequence

<img style="display: inline-block" src="/screenshots/raccoon_main.png" width="400" height="400">


| Parameter | Default       |   Description |
| :---      | :---          | :---          |
| -s        | seq.txt       | Filename with extension for the sequence.                            |
| -o        | out           | Filename without extension for the resulting pdb file.               |
| -m        | monomers.dat  | Filename with extension for the monomer library.                     |
| -e        | False         | Write all explicit bonds. It can be useful for nonstandard residues. |
| -r        | True          | Remove duplicate bonds. Bond a --> b is equivalent to b --> a.       |

### Check PDB Files & Export XYZ Files

## Limitations

> [!WARNING]
> * The generated output is only a starting structure for molecular dynamics simulations.
> * Users must ensure atom/bead naming consistency for matching force field data.

