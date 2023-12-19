###### Functions
import numpy as np
import os
import pandas as pd
import ast
import copy
# Initialize monomers from a text file and import them to the python environment.

def InitializeMonomers(file_name):
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

#  Normalization of Atom Positions

def Normalization(monomer):

    reference = [0] * 3
    for x in range(0,3):
        reference[x] = monomer[1][x+2]
    for x in range(1,(monomer['atoms']+1)):
        for y in range(2,5):
            monomer[x][y] = np.round((monomer[x][y] - reference[y-2]), 5)
    return (monomer)


# Updating the monomer position after each succesful monomer addition step.

def UpdateMonomer(monomer, shift, shift_cartesian):
    
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

    
# This function opens the PDB file and adds a polymer chain with a given repeating unit, chain length and resid_shift
# Might be a useful function for long repeating chains

def CreatePolymerChain(monomer, chain_length, filename, resid_shift):
    
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

# Main function of the script
# adds monomers to a polymer peptide chain
def AddMonomer(monomer, filename, resid, shift):
    
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

# This function is needed to sort the PDB file after adding the bonds               
def Sorting(filename):
    with open(f"{filename}.pdb", "r") as input_file, open(f"{filename}-copy.pdb", "w") as output_file:
        rows = input_file.readlines()
        sorted_rows = sorted(rows, key=lambda x: x.split()[0])
        output_file.writelines(sorted_rows)

# Closes the PDB file. Needs the overall number of atoms / beads!
def CloseFile(filename, number):
     with open(f"{filename}.pdb", "a") as file:
        file.write("{:>0}{:<11}{:<5}{:<5}{:<4}{:<5}{:<5}{:<5}{:<5}{:<5}{:<5}{:<5}{:<5}{}\n".format("", "MASTER", 0, 0, 0, 0, 0, 0, 0, 0, number, 0, number, 0))
        file.write("END")

# Creates the lines for the inter-monomer bonds / linkages. Needs a list of all links.
def Bonds(links, filename):
    for i in range(0, len(links)-1, 1):
        with open(f"{filename}.pdb", "a") as file:
            file.write("{:>0}{:<7}{:<5}{:<5}{:<4}{:<3}{:<6}{:<8}{:<8}{:<10}{:<7}{:<14}{}\n".format("", "CONECT", links[i], links[i+1], "", "", "", "", "", "", "", "", "", ""))

def Sequence(file_name):
    
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

# Generates the pdb file based on a sequence (text file input)
# Inputs: sequence file (e.g. sequence.txt), monomer file (e.g. monomers.dat), output file (e.g. my_pdb_file.pdb)

def GenerateFile(sequence_file='seq.txt', monomer_file='monomers.dat', out_file='outfile', explicitbonds = False, remove_duplicates=True, debug=True):
    
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
    return copy.deepcopy(x)

def InvertAminoAcid(monomer):
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
