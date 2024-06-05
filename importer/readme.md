## Simplified PDB Import

Using the `importer.py` script lets you add PDB structures directly to the `monomers.dat` file. Please run the command as follows:
```bash
python importer.py -f yourpdb.pdb >> ../monomers.dat
```
After adding the monomers to the `monomers.dat` file please use the following two commands in a Jupyter notebook to update the available monomer library:
```python
import project_raccoon as rc 
monomers = rc.Monomers.from_file('monomers.dat')
```
We implemented the possibility of atom name mapping so that the monomer import for different force fields becomes easier. Our provided *opls.map* file maps the standard [Avogadro](https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-4-17) PDB naming scheme to the [OPLS-AA](https://pubs.acs.org/doi/10.1021/acs.jctc.5b00356) naming scheme for amino acids. The script supports the following input parameters:

| Parameter | Default       |   Description |
| :---      | :---          | :---          |
| -f        | -             | Filename with PDB extension for the import. |
| -m        | opls.map      | Filename with extension for atom name mapping file.|
| -c        | 'C'           | Atom name for the C terminus.|
| -n        | 'N'           | Atom name for the N terminus. |
| -p        | 0             | Saves the monomer as polymer (1) or not (0).|
| -r        | 'atomistic'   | Resolution of the monomer.|

Please remember you add new monomers to your *residuetypes.dat* file and provide an appropiate topology file for *GROMACS*. Always ensure the right atom names. Sometimes the `pdb2gmx -glu -his` etc. flags are necessary to provide the right protonation state.

> [!NOTE]
> **Transparency Statement**: We declare that the code sections for the import function were written with [Open AI GPT-4](https://doi.org/10.48550/arXiv.2303.08774) assistance.
