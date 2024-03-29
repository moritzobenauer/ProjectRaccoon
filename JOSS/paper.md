---
title: 'Project RACCOON: Automated construction of PDB files for polymers and polymer peptide conjugates'
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
    equal-contrib: true
  - name: Kai N. Spauszus
    orcid: 0009-0006-0650-2273
    affiliation: 1
    equal-contrib: true
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

The *Project RACCOON* software was developed to automatically generate PDB files for linear polymers, polypeptides, and polymer peptide conjugates. Previously published software cannot easily represent polypeptide strands with converging termini. *Project RACCOON* also offers the advantage of combining multiple modeling resolutions within one system. Although this software was explicitly developed for C2 symmetric polymer peptide conjugates introduced by @otter2018 with the OPLS-AA/M force field, its application is not limited to this molecule class. [@oplsaa2019; @oplsaa_proteins] PDB files can be created for all polymers, polypeptides with natural and non-natural amino acids, and other macromolecules combined with any force field.

# Statement of need

Common software packages for generating the input structures of molecular dynamics simulations such as *polyply* or *CHARMM-GUI* are unsuitable for representing the class of C2 symmetric polymer peptide conjugates. [@polyply; @charmm-gui] This class of molecules is characterized by the convergence of two C-termini in the polypeptide sequences. This represents a non-natural system and has not yet been implemented in the aforementioned software packages. Modeling these sequences is straightforward and achievable for novice users with *Project RACCOON*, as the flipping of N- and C-termini of amino acids is an integral part of the software.
Construction of the molecules is carried out in a building block fashion, allowing the automatic import of natural and non-natural amino acid repeating units and synthetic polymers. In addition, this modular approach offers the advantage of combining building blocks with different resolutions (atomistic, united atom, and coarse-grained) in one molecule. This is necessary to describe polymer peptide conjugates on hybrid scales. [@taylor2020]
The software presented here is an improved alternative to other options like *polyply* because the generated files are directly available in the PDB standard format, which can be processed and visualized by any cheminformatic software. Until now, the PDB format has been used almost exclusively for proteins in the bioinformatics context. [@pdb-files] However, this work shows that extending the PDB files to general macromolecules simplifies and standardizes data processing. This software represents an important step towards standardization and simplification in the field of cheminformatics and thus contributes to the idea of the UNESCO Open Science Initiative. [@open-science]

# Functionality and Extensibility

New monomers can be quickly and systematically added to the software. Each monomer is assigned various properties such as the number of atoms, linkage points, and the geometry of the monomer building block. This data is saved as a JSON file and can therefore be easily imported, edited, and exported. [@python3] The geometry of the monomer building block can be obtained from vacuum energy-minimized structures or crystallographic data. Bonds link all monomers together, which the user can also define, and are only limited by the capacity of the force field used later to account for these bonds. Each atom is also assigned an individual name, which serves as an identifier for the associated force field.
Novice users can utilize the automated import feature on molecular building blocks generated with software like *Avogadro 2*. [@avogadro] Advanced users can import *Project RACCOON* into their projects and use all functions for further future work, as the software can quickly adapt to their needs. Functionality and stability were extensively verified with unit testing, resulting in a code coverage of over 85%. *Project RACCOON* can enhance productivity and collaboration for all users since one can build libraries of monomers and share them with others.
It should be noted, however, that the geometry is generated from a self-avoiding random walk and has no physical justification. Parameters for this self-avoiding random walk can be chosen by the user to obtain a pseudo phantom chain structure. Energy minimization steps and unrestricted equilibration steps must be performed with suitable software in any case before using the model in molecular dynamics simulations.

# Acknowledgements

Funding from the DFG (*Deutsche Forschungsgemeinschaft*) is acknowledged. P.B. and F.S. are members of the GRK 2516 (Project No. 405552959). This work was supported by the German Academic Exchange Service (*DAAD*) International Study and Training Partnerships (*ISAP*) program (Grant No. 57460493). M.L.O. received a scholarship to pursue a research internship in the research group of M. Muthukumar, Polymer Science and Engineering, University of Massachusetts, Amherst. 

# Conflict of Interest

The authors declare no conflict of interest.

# References
