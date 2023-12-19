# Project RACCOON
**Automated construction of atomistic and coarse-grained models in the PDB format for linear polymer peptide conjugates.**

## General Purpose & Scope
 <div align="justify">Project RACCOON (*Rapid Automated Construction of Conjugates using the Optimized OPLS Input*) is a Python tool designed for the generation of PDB (Protein Data Bank) files for polymer peptide conjugates, polypeptides, and polymers in a building block fashion. It allows for the easy addition of new monomers, incorporation of polypeptide and polymer sequences in text form, and outputs PDB and XYZ files.
The tool was specifically developed for the OPLS force field but is adaptable to any other force field. The tool's primary scope is to provide a straightforward method for creating starting structures for molecular dynamics simulations. It caters to the creation of linear polymer peptide conjugates with facing polypeptide strands. 
The tool emphasizes user-friendliness, making it accessible to chemists and physicists, not just limited to theoretical chemists with extensive background knowledge. Additionally, users can modify the code for their specific projects. Early developments started in January 2023. The software was developed at the KOMET and Besenius' research groups at the Johannes Gutenberg University Mainz.</div>

## How to use RACCOON

### Importing new Monomers

<img src="/screenshots/import1.png" width="100" height="100">

## Limitations

> [!WARNING]
> * The generated output is only a starting structure for molecular dynamics simulations.
> * Users must ensure atom/bead naming consistency for matching force field data.

