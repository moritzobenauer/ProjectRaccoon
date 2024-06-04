## Geometry Builder

> [!CAUTION]
> This code is not peer-reviewed yet and is considered to be error-prone. Since it is necessary for current research projects, we decided to share it anyway as part of the **Project RACCOON** package. 

After generating different macromolecules with **Project RACCOON** it can be necessary to align and combine different molecules e.g. to simulate protein-ligand systems. We provide an interactive geometry builder that takes *pdb* files as an input, visualizes them using [**py3Dmol**](https://github.com/avirshup/py3dmol) and lets the user translate and rotate the molecules. After the desired geometry is achieved, the updated *pdb* files can be saved.  

During this development stage the main **Project RACCOON* and *Project RACCOON* Geometry Builder dependencies are not harmonized yet. Therefore, we recommend to use a seperate python environment to run *geometry_builder.py*. Installing the necessary dependencies using [mamba](https://github.com/mamba-org/mamba):

```bash
mamba create --name prdev --file requirements.txt
mamba activate --name prdev
```
After the succesful installation you can run the program with the following command: 

```bash
python3 geometry_builder.py -f file1.pdb file2.pdb
```

The following keyboard shortcuts are used in the standard configuration.

Keyboard Character | Function
-------- | -------- 
<kbd>Q</kbd> <kbd>E</kbd>    | Flip Molecule    
<kbd>W</kbd> <kbd>A</kbd> <kbd>S</kbd> <kbd>D</kbd>    | Rotate Molecule
<kbd>U</kbd> <kbd>N</kbd>    | Move Molecule in up and down direction
<kbd>H</kbd> <kbd>J</kbd>    | Move Molecule in right and left direction
<kbd>C</kbd>                 | Change currently selected molecule
<kbd>1</kbd>                 | Close without saving new coordinates
<kbd>2</kbd>                 | Close with saving new coordinates


**Project RACCOON** Geometry Builder uses the <kbd>F5</kbd> key to reload web pages. If your standard browser does not support that please change the source code accordingly. The standard rotation angle $\Theta$ and standard translational increment *dx* can be set with the `-t` and `-m` flag respectively. Please use the `--help` flag for more information.


> [!NOTE]
> **Transparency Statement**: We declare that the code sections for the Geometry Builder were written with [Open AI GPT-4](https://doi.org/10.48550/arXiv.2303.08774) assistance.
