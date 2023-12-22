---
title: 'Project RACCOON: Extending the PDB file format to polymers and polymer peptide conjugates'
tags:
  - Python
  - cheminformatics
  - molecular-dynamics-simulatiom
  - pdb-files
  - molecular-modeling
authors:
  - name: Moritz L. Obenauer
    orcid: 0009-0003-8140-9907
    affiliation: 1 
    corresponding: true 
  - name: Kai N. Spauszus
    orcid: 0009-0006-0650-2273
    affiliation: 1 
- name: Pol Besenius
    orcid: 0000-0001-7478-4459
    affiliation: 1 
- name: Friederike Schmid
    orcid: 0000-0002-5536-6718
    affiliation: 2 
affiliations:
 - name: Department of Chemistry, Johannes Gutenberg-University Mainz, Duesbergweg 10-14, D-55128 Mainz, Germany
   index: 1
 - name: Institute of Physics, Johannes Gutenberg University, Staudingerweg 9, 55128 Mainz, Germany
   index: 2
date: 20 January 2024
bibliography: paper.bib
---

# Summary

The Project RACCOON software was developed to automatically generate PDB files for linear polymers, polypeptides and polymer peptide conjugates. Previous published software cannot easily represent two converging termini of polypeptide strands. Besides this unique feature, Project RACCOON offers the advantage of combining different modelling resolutions in one system. Although this software was explicitly developed for C2 symmetric polymer peptide conjugates introduced by Otter et. al. with the OPLS-AA/M force field, the application is not limited to it. [@oplsaa2019, @oplsaa_proteins, @otter2018] PDB files can be created for all types of polymers, polypeptides with natural and non-natural amino acids, and other macromolecules in combination with any force field.

# Statement of need

The basis for developing this software was the problem that classic software packages such as polyply or CHARMM-gui for generating input structures for molecular dynamics simulations are not suitable for representing the class of C2 symmetric polymer peptide conjugates.[@polyply, @charmm-gui] These are characterized by the two C-termini of the polypeptide sequences converging. This represents a non-natural system and has not yet been implemented in these software packages. Modeling these sequences is very easy with Project RACCOON, as the inversion of amino acids is an integral part of the software.
Modeling is carried out in a building block fashion. This allows the automatic import of synthetic polymers and natural and non-natural amino acid repeating units. In addition, this offers the advantage of combining building blocks with different resolutions (atomistic, united atom, coarse-grained) in one molecule. This is necessary in order to be able to describe polymer peptide conjugates on hybrid scales.[@taylor2020]
The software presented here is an alternative to other options like polyply because the generated files are available in the PDB standard format, which can be processed and visualized by any software. Until now, the PDB format has been used almost exclusively for proteins in the bioinformatics context.[@pdb-files] However, this work shows that extending the PDB files to general macromolecules simplifies and standardizes data processing. This software represents an important step towards standardization and simplification in the field of cheminformatics and thus contributes to the idea of the UNESCO Open Science Initiative.[@open-science]

# Extensibility

New monomers can be quickly and systematically added to the software. Each monomer is assigned various properties in a Python dictionary, such as the number of atoms, linkage points, and the geometry of the monomer building block.[@python3] The geometry of the monomer building block can be obtained from vacuum energy-minimized structures or crystallographic data. Bonds link all monomers together, which the user can also define, and are only limited by the capacity of the used force field to account for these bonds. Each atom is also assigned an individual name, which serves as an identifier for the force field.
Novice users can use the automated import of molecular building blocks generated with software like Avogadro 2.[@avogadro] Advanced users can use the source code as a reference or starting point for further projects, as the software can quickly adapt to their needs. Furthermore, one can build libraries of monomers and share them with other users.
It should be noted, however, that the generated structure has no physical justification, and energy minimization steps and unrestricted equilibration steps must be performed in any case before using the model in molecular dynamics simulations.

# Acknowledgements

The authors acknowledge funding through the Deutsche Forschungsgemeinschaft (DFG) through GRK 2516 (Grant No. 405552959). This work was supported by the German Academic Exchange Service. M.L.O. was the recipient of a scholarship to pursue a research internship in the research group of M. Muthukumar, Polymer Science and Engineering, University of Massachusetts, Amherst. 

# Conflict of Interest

The authors declare no conflict of interest.

# References
