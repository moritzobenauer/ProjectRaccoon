"""
Beschreibung worum es geht. lizenz etc...

"""

import numpy as np              # erkläre wofür imports genutzt werden
import os
import pandas as pd
import ast
import copy
from functions_advanced import *
from biopandas.pdb import PandasPdb
ppdb = PandasPdb().fetch_pdb('3eiy')
from itertools import combinations
import argparse
from rich.console import Console
from questionary import *

argParser = argparse.ArgumentParser()
argParser.add_argument("-s", "--sequence", help="Sequence File", type=str, default="seq.txt")
argParser.add_argument("-o", "--output", help="Output File", type=str, default="out")
argParser.add_argument("-m", "--monomers", help="Monomer File", type=str, default="monomers.dat")
argParser.add_argument("-e", "--explicit", help="Explicit Bonds", type=bool, default=False)
argParser.add_argument("-r", "--removeduplicates", help="Remove Duplicates", type=bool, default=True)

args = argParser.parse_args()

sequencefile = args.sequence
outputfile = args.output
monomerfile = args.monomers
explicitbondsbool = args.explicit
removeduplicatesbool = args.removeduplicates


def InitializeMonomers(file_name:str):

    """Reads inputsfile and parses into results?
        Frage: Was ist results? Wie sieht der Input aus? Wie ist die Struktur von Results?_> Ducktyping

    Args:
        file_name (str): _description_

    Returns:
        _type_: _description_
    """
    result = [{}]
    with open(file_name) as inpt:
        for line in inpt:
            if line.startswith('#'):
                continue
            if line.strip() == '':
                result.append({})
            else:
                key, value = line.split('=')
                result[-1][key] = value.strip()
        for item in result:
            item['polymer'] = bool(int(item['polymer']))
            item['atoms'] = int(item['atoms'])
            item['link'] = ast.literal_eval(item['link'])
            for x in range(1,item['atoms']+1):
                item[x] = ast.literal_eval(item.pop(f'{x}'))
            item['inverted'] = False
        amino_acids = [d for d in result if d.get('polymer', False) is False]
        for item in amino_acids:
            # Some capping groups cannot be inverted
            try:
                result.append(InvertAminoAcid(item))
            except:
                pass
    return result

def Normalization(monomer):

    """Normalization of Atom Positions
    Frage: Normalisierung aller Atome im Molekül??

    Args:
        monomer (_type_): _description_

    Returns:
        _type_: _description_
    """

    reference = [0] * 3
    for x in range(0,3):
        reference[x] = monomer[1][x+2]
    for x in range(1,(monomer['atoms']+1)):
        for y in range(2,5):
            monomer[x][y] = np.round((monomer[x][y] - reference[y-2]), 5)
    return (monomer)


def UpdateMonomer(monomer, shift, shift_cartesian):

    """Updating the monomer position after each succesful monomer addition step.

    Args:
        monomer (_type_): _description_
        shift (_type_): _description_
        shift_cartesian (_type_): _description_
    """
    
    # Update the linking atoms / beads due to the renaming of the atoms / beads
    
    try:
        monomer["link"][0] = monomer["link"][0] + (shift-1)
        monomer["link"][1] = monomer["link"][1] + (shift-1)
    except:
        pass
    
    i = 1
    
    while i <= monomer["atoms"]:
        
        # Update Cartesian Positions
        
        for j in range(2,5):
            monomer[i][j] = np.round(monomer[i][j] + shift_cartesian[j-3], 2)
            
        # Renaming the atoms / beads
        
        monomer[i][6] += (shift-1)
        
        # Renaming the neighboring atoms / beads
        
        for j in range(0, len(monomer[i][5])):
            monomer[i][5][j] = monomer[i][5][j] + (shift-1)
        i += 1  

    
def CreatePolymerChain(monomer, chain_length, filename, resid_shift):

    """This function opens the PDB file and adds a polymer chain with a given repeating unit, chain length and resid_shift
        Might be a useful function for long repeating chains.

    Args:
        monomer (_type_): _description_
        chain_length (_type_): _description_
        filename (_type_): _description_
        resid_shift (_type_): _description_

    Return:
        noting
    """
    
    global links
    global atom_counter
    
    links = []
    atom_counter = monomer['atoms'] * chain_length
    
    # Ensures a conversion from standard cartesian coordinates to coordinates that are centered around the first atom / bead
    
    monomer = Normalization(monomer)
    
    for unit in range(chain_length):
        # Setting the resid (this can be shifted)
        resid = unit + 1 + resid_shift
        
        # Adding linkage points of the monomer to a global list
        
        # Tries to add the C-Terminus links
        try:
          links.append(monomer["link"][1])
        except:
          pass
        # Tries to add the N-Terminus links
        try:
          links.append(monomer["link"][0])
        except:
          pass
        
        with open(f"{filename}.pdb", "a") as file:
            i=1
            # Writing every atom / bead of the monomer in a new line
            while i <= monomer["atoms"]:
                # Ensures right PDB file construction (spaces are important!)
                file.write("{:>0}{:<7}{:<5}{:<5}{:<4}{:<3}{:<6}{:<8}{:<8}{:<10}{:<7}{:<14}{}\n".format("", "ATOM", monomer[i][6], monomer[i][0], monomer["res"], "A", resid, monomer[i][2], monomer[i][3], monomer[i][4], 1.0, 0.0, monomer[i][1]))
                i += 1
        # UpdateDic: Updating the monomer for the next chain segment
        # s: Renames the atoms / beads consecutively
        # m: list of xyz cartesian shifts (this could be a function itself, so it is not only linear)
        m = [1,1,1]
        s = monomer['atoms'] + 1
        UpdateMonomer(monomer, s, m)
        file.close()      

def AddMonomer(monomer, filename, resid, shift):

    """
    Central function of the modul: adds monomers to a polymer peptide chain.

    Args:
        monomer (_type_): _description_
        filename (_type_): _description_
        resid (_type_): _description_
        shift (_type_): _description_
    """
    
    global links
    global atom_counter
    global res_counter
        
    atom_counter += monomer['atoms']
    res_counter += 1
    
    # Ensures a conversion from standard cartesian coordinates to coordinates that are centered around the first atom / bead
    
    monomer = Normalization(monomer)
    
    # Shifting cartesian coordinates and linking points to the right position
    
    # UpdateMonomer: Updating the monomer for the next chain segment
    # s: Renames the atoms / beads consecutively
    
    # m: list of xyz cartesian shifts
    # The r_min, r_max values can be changes with their respective functions. These should only be changed in case of an overlap of two atoms which can lead to infinite forces.
    # The output structure is never physically correct and need to be energy minimized in all cases.
    # chain is elongated along z direction

    global x_min
    global x_max
    global y_min
    global y_max
    global z_min
    global z_max

    x_min, x_max = -5.5, 5.5
    y_min, y_max = -5.5, 5.5
    z_min, z_max =  0.3, 0.4
    
    m = np.array([(np.random.uniform(x_min, x_max)),(np.random.uniform(y_min, y_max)),(shift * np.random.uniform(z_min, z_max))])

    s = shift + 1
    UpdateMonomer(monomer, s, m)
    
    
    # Adding linkage points of the monomer to a global list
        
    # Tries to add the C-Terminus links
    # TODO: so wenig try expect wie möglich. Weshalb ist das hier nötig? Anderweitig abfangen 
    try:
        links.append(monomer["link"][1])
    except:
        pass
    # Tries to add the N-Terminus links
    try:
        links.append(monomer["link"][0])
    except:
        pass
        
    with open(f"{filename}.pdb", "a") as file:
        i=1
        # Writing every atom / bead of the monomer in a new line
        while i <= monomer["atoms"]:
            # Ensures right PDB file construction (spaces are important!)
            file.write("{:>0}{:<7}{:<5}{:<5}{:<4}{:<3}{:<6}{:<8}{:<8}{:<10}{:<7}{:<14}{}\n".format("", "ATOM", monomer[i][6], monomer[i][0], monomer["res"], "A", resid+1, monomer[i][2], monomer[i][3], monomer[i][4], 1.0, 0.0, monomer[i][1]))
            i += 1
        file.close()      

# Easy access to the r_min, r_max variables
# TODO: Unschöner Stil: eher globale variablen definieren, müssen die global sein? reicht nicht beispielsweise auch ein struct, wie named tuple?

def SetXMinMax(min, max):
    
    global x_min
    global x_max
    x_min, x_max = min, max

def SetYMinMax(min, max):
    global y_min
    global y_max
    y_min, y_max = min, max

def SetZMinMax(min, max):
    global z_min
    global z_max
    z_min, z_max = min, max
              
def Sorting(filename:str):
    """This function is needed to sort the PDB file after adding the bonds 

    Args:
        filename (str): _description_

    Returns:
        nothing
    """

    

    with open(f"{filename}.pdb", "r") as input_file, open(f"{filename}-copy.pdb", "w") as output_file:
        rows = input_file.readlines()
        sorted_rows = sorted(rows, key=lambda x: x.split()[0])
        output_file.writelines(sorted_rows)

def CloseFile(filename, number):
    """Closes the PDB file. Needs the overall number of atoms / beads!
    # TODO: format schöner schreiben

    Args:
        filename (_type_): _description_
        number (_type_): _description_
    """
    with open(f"{filename}.pdb", "a") as file:
        file.write("{:>0}{:<11}{:<5}{:<5}{:<4}{:<5}{:<5}{:<5}{:<5}{:<5}{:<5}{:<5}{:<5}{}\n".format("", "MASTER", 0, 0, 0, 0, 0, 0, 0, 0, number, 0, number, 0))
        file.write("END")

def Bonds(links, filename):
    """Creates the lines for the inter-monomer bonds / linkages. Needs a list of all links.

    Args:
        links (_type_): _description_
        filename (_type_): _description_

    Returns:
        nothing
    """
    for i in range(0, len(links)-1, 1):
        with open(f"{filename}.pdb", "a") as file:
            file.write("{:>0}{:<7}{:<5}{:<5}{:<4}{:<3}{:<6}{:<8}{:<8}{:<10}{:<7}{:<14}{}\n".format("", "CONECT", links[i], links[i+1], "", "", "", "", "", "", "", "", "", ""))

def Sequence(file_name:str):
    """_summary_

    Args:
        file_name (str): _description_

    Returns:
        _type_: _description_
    """
    
    result = []
    sequence = []
    with open(file_name) as inpt:
        for line in inpt:
            if line.startswith('#'):
                continue
            if line.strip() == '':
                result.append({})
            else:
                res, resolution, inverse, repeat = line.split(':')
            
                
                # Look-Up Tabble for resolution
                for x in range(int(repeat)):
                    if resolution == 'AA':
                        resolution_lookup = 'atomistic'
                    elif resolution == 'UA':
                        resolution_lookup = 'united_atom'
                    elif resolution == 'CG':
                        resolution_lookup = 'coarse_grained'
                    index = next((i for i, d in enumerate(monomers) if (d.get('res') == res) and (d.get('resolution') == resolution_lookup) and (d.get('inverted') == bool(int(inverse)))), None)
                
                    result.append(f'monomers[{index}]')        

    return result

def GenerateFile(sequence_file='seq.txt', monomer_file='monomers.dat', out_file='outfile', explicitbonds = False, remove_duplicates=True, debug=True):
    """Generates the pdb file based on a sequence (text file input)
       Inputs: sequence file (e.g. sequence.txt), monomer file (e.g. monomers.dat), output file (e.g. my_pdb_file.pdb)

       # TODO: globale Variablen 

    Args:
        sequence_file (str, optional): _description_. Defaults to 'seq.txt'.
        monomer_file (str, optional): _description_. Defaults to 'monomers.dat'.
        out_file (str, optional): _description_. Defaults to 'outfile'.
        explicitbonds (bool, optional): _description_. Defaults to False.
        remove_duplicates (bool, optional): _description_. Defaults to True.
        debug (bool, optional): _description_. Defaults to True.
    """
        
    global debugging

    # Prevent data loss

    try:
        os.rename(f'{out_file}.pdb',f'{out_file}_old.pdb')
    except:
        pass
    
    global atom_counter
    global res_counter
    global links
    global monomer_bonds
    global filtered_list
    global removeduplicates
    global debugging
    
    
    atom_counter, res_counter = 0,0
    monomer_bonds,filtered_list, links = [],[], []
    
    # The bonds (x,y) are equal to (y,x)
    
    removeduplicates = remove_duplicates
    
    
    
    for item in Sequence(sequence_file):
        monomers = InitializeMonomers(monomer_file)
        if explicitbonds == True:
            ExplicitBonds(eval(item))
        else:
            pass
        AddMonomer(eval(item), out_file, res_counter, atom_counter)
        
    # After iterating through all monomers bonds are written to the file, lines are getting sorted and file gets closed with PDB closing statement   
    
    WriteExplicitBonds(out_file)
    Bonds(links, out_file)
    Sorting(out_file)
    CloseFile(f'{out_file}-copy', atom_counter)
    os.remove(f'{out_file}.pdb')
    os.rename(f'{out_file}-copy.pdb',f'{out_file}.pdb')
    #

    print(f'Created PDB file with {res_counter} residue, {atom_counter} atoms, and {int(len(links))} inter-residual bonds.')
    
    # Resets standard values
    
    atom_counter = 0
    res_counter = 0
    links = []


def GetCopy(x):
    # TODO: Brauchst du hierfür eine Funktion?
    return copy.deepcopy(x)

def InvertAminoAcid(monomer):
    """_summary_

    Args:
        monomer (_type_): _description_

    Returns:
        _type_: _description_
    """
    monomer_inverted = GetCopy(monomer)
    temp = GetCopy(monomer_inverted['link'])
    monomer_inverted['link'][0] = temp[1]
    monomer_inverted['link'][1] = temp[0]
    monomer_inverted['inverted'] = True
    return monomer_inverted


def ExplicitBonds(monomer):
    
    global filtered_list
    
    seen_pairs = set()
    for x in range(atom_counter,atom_counter+monomer['atoms']):
        for y in range(len(monomer[x-atom_counter+1][5])):
            monomer_bonds.append((x+1, monomer[x-atom_counter+1][5][y]+atom_counter))
            
    if removeduplicates == True:
        # This function eliminates duplicates
        seen_pairs = set()
        filtered_list = []       
        for pair in monomer_bonds:
            reversed_pair = (pair[1], pair[0])
            if reversed_pair not in seen_pairs:
                seen_pairs.add(pair)
                filtered_list.append(pair)
    else:
        filtered_list = monomer_bonds

        
def WriteExplicitBonds(out_file):
    # Writing to the outfile
    for item in filtered_list:
        with open(f"{out_file}.pdb", "a") as file:
            file.write("{:>0}{:<7}{:<5}{:<5}{:<4}{:<3}{:<6}{:<8}{:<8}{:<10}{:<7}{:<14}{}\n".format("", "CONECT", item[0], item[1], "", "", "", "", "", "", "", "", "", ""))



def Welcome():
    # TODO: auslagern?
    console.print("""                                                                                                     
                                                                                                    
            ..            ..             ..            ..             ..            ..             .
     *******//*************/*************//*************/*************/*************//*************/
                                                                                                    
            ..            ..             ..            ..             ..            ..             .
     ,,,,,,,**,,,,,,,,,,,,**,,,,,,,,,,,,,**,,,,,,,,,,,,**,,,,,,,,,,,,,**,,,,,,,,,,,,**,,,,,,,,,,,,,*
                                                                                                    
            ..            ..             ..        /   ..    ,#(.     ..            ..             .
     .......,,............,,,.........,*/(((##/*,.*/((*,**/#(/(/.....,,,............,,............,,
                                     /&&%%#(((/((/*,,,,,,,,.,*#@                                    
            ..            ..        .##//(#%%&&%#(*****///*******#(   ..            ..             .
           ...            ...    .%(//*//(///////#%#//////*////*///((...       ,*.  ...           ..
                             ,&%..,,.,,,,*/#%(,,.,,,,,.             ,(        ,     ./*             
            ..            .(#/*********//((/******,    ..   .*%&/&&&&(..      /..**//*, ./.        .
            ..          #%%////*/*//*#%(/////*.       ,(&&&&&/   %&&&%(*     .,..***//(///.,/      .
                           ,*,,*,,,.*#/     ./@@@@&&&&&&%%%%&%&&#,      (#.  *  .......,... ,       
            ..           *//************//(#%&&&&&&&&&&&&&&&#,   //*,*,.     *..*,,,***,,,*,*      .
            ..          ((///*//**//*///////***//(**/*/(#(#&(((##%%/. ..    ,..,**,****,,,*,/      .
                       &*/***,.,,/*..,,,*/***.*,,,*,*,******/%              /  ,......,.*.. ,       
            ..        %%(**///**********///**/*/****/**///**/(/       ..   .,..*,,,,****(**,*      .
            ..       #(*//////*//***/**///////***/*////*///**//%      ..   *..,*,,*,***/*,*,*      .
                   ,%%@&***/**,,.,/,.   ,*(#@@@.            .*** (@        /  *.....,,,..., ,       
            ..        %***///*********/*///*****%(     *((/. ****///##*,   #...,/(/****,,,*,*      .
            ..       ((/**/(//***/***/*////////**/%*   (,   .#%((/*#&**,   .#(*,....,,.,,//.,      .
            .       (/&%/,****,**,/,,,****//,***/*,*@  (,                             ,* ,%         
            ..        .#/**/************///********/%  /,.            ..           .*,,,           .
            ..         &#////****//**/**///***//*/*/&  ./*     ,***,,,***/*.       .*,,,           .
     ,......,,,,........,@#(/*(*/****(/*////(*//*,#& ..,((*%,*,,,*,,*,,*************,,,,*******(..,,
            ..           ../(************//*****/%%(/****************(((###/********************   .
            ..           /#/////(///(//(((##%&#(*******//************(##(((*********///********/   .
     *****,*********,*@@(*/////#.**,,,,,****,******,,*,***,,*,**,*,******,**,***,*,****,************
            ..      #/%&&&(/*((          ..            ..             ..            ..             .
            ..    (#*/**%&&&%(           ..            ..             ..            ..             .
     /////*////#@ &%@@((/##  *////////////////////////////////////////////////////////////*/////////
            .*/*/**(&&##. ..             ..            ..             ..            ..             .
           .(*,*###/.      .             ..            ..             ..            ..             .
     (((/(/(((((((/(///(/((((/(((/((((//((((((((/((((//(((/(//(/((((/((((/(/(((((((/(((/((//((((//((
            ..            ..             ..            ..             ..            ..             .
            ..            ..             ..            ..             ..            ..             .
""")
    
    console.print('Project RACCOON by Moritz Obenauer @ JGU Mainz 2023 \n \n')
    
def main():
    global monomers 
    
    options = select("Choose Function",choices=["Create PDB File", "Check PDB File", "Convert PDB to XYZ File","Check Minimal Distance","Exit"],).ask()

    while True:

        if options == 'Create PDB File':
            monomers = InitializeMonomers('monomers.dat')
            GenerateFile(sequence_file=sequencefile, monomer_file=monomerfile, out_file=outputfile, explicitbonds = explicitbondsbool, remove_duplicates=removeduplicatesbool)
            main()
        elif options == 'Check PDB File':
            CheckPDB(outputfile)
            main()
        elif options == "Convert PDB to XYZ File":
            PDBtoXYZ(outputfile)
            main()
        elif options == 'Check Minimal Distance':
            CheckMinimalDistance(outputfile)
            main()
        elif options == 'Exit':
            exit()

if __name__ == "__main__":
    console = Console()
    Welcome()
    main()
