# Import Avogadro 2 .bs files as new monomers

# Setting up environment

import pandas as pd
import argparse
from tabulate import tabulate
pd.set_option('display.width', None)
tabulate.PRESERVE_WHITESPACE = True
from rich.console import Console
from questionary import *
import os

def clear_terminal():
    os.system('cls' if os.name == 'nt' else 'clear')

# User can set a bs input file and monomer output file

argParser = argparse.ArgumentParser()
argParser.add_argument("-f", "--input", help="Input File", type=str)
argParser.add_argument("-o", "--output", help="Output File", type=str)

args = argParser.parse_args()

inputfile = args.input
outputfile = args.output

# Function to write to the monomers.dat file in the right format

def AddToDat():
    new_lines = ['\n',
                 'res=' + structure_name + '\n', 
                 'resolution=' + resolution + '\n',
                 'polymer=' + str(polymer) + '\n',
                 'atoms=' + str(len(structure)) + '\n',
                 'links=' + str(links) + '\n']
    
    with open(outputfile, 'a') as file:

        file.writelines(new_lines)

        for index, row in structure.iterrows():
            line = f'''{row[0]}=["{row[1]}", "{row[2]}", {row[4]}, {row[5]}, {row[6]}, {row[3]}, {row[0]}] \n'''
            file.write(line)

        file.close()

    console.print(f'Monomer added to {outputfile}')

def main():

    global structure_name
    global polymer
    global resolution
    global links
    global structure

    header = ["Element", "x", "y", "z", "N1", "N2", "N3", "N4", "#", "Neighbors", "Identifier"]

    # Import the bs file, sep is four spaces

    structure = pd.read_csv(inputfile, sep='   ', skiprows=[0,1], header=None, names=header)

    def combine_columns(row):
        values = [int(value) for value in [row['N1'], row['N2'], row['N3'], row['N4']] if not pd.isna(value)]
        return values
    def update_column_hash(row):
        return row.name + 1  # Adding 1 to the row index
        
    structure['Neighbors'] = structure.apply(combine_columns, axis=1)
    structure['#'] = structure.apply(update_column_hash, axis=1)
    structure = structure[["#", "Identifier", "Element", "Neighbors", "x", "y", "z"]]

    clear_terminal()

    structure_name = text(f"Enter a name for the new monomer").ask()

    polymerQ = confirm("Is this a standard polymer?").ask()
    if polymerQ == True:
        polymer = 1
    else:
        polymer = 0

    resolution = select("Choose resolution",choices=["atomistic", "united_atom", "coarse_grained"],).ask()

    for index, row in structure.iterrows():
        user_input = text(f"Enter a force field Identifier for {row['Element']}' with the Number {row['#']}: ").ask()
        structure.at[index, 'Identifier'] = user_input
    clear_terminal()
    console.print(tabulate(structure))

    linkC = text(f"Choose C-Terminus (1-{len(structure)})").ask()
    linkN = text(f"Choose N-Terminus (1-{len(structure)})").ask()

    links = [int(linkC), int(linkN)]

    save = confirm("Save the monomer to the monomer data file?").ask()

    if save == True:
        AddToDat()
    else:
        pass

# code is written within a main function for interactive user input

if __name__ == "__main__":
    console = Console()
    main()

