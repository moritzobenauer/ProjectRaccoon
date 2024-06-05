import argparse
import pandas as pd

def parse_pdb(file_path, atom_mapping):
    atoms = []
    conect_records = []
    global atom_info

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("ATOM"):
                atom_name = line[12:16].strip()
                if atom_name in atom_mapping:
                    atom_name = atom_mapping[atom_name]
                atom_info = {
                    "atom_serial": int(line[6:11].strip()),
                    "atom_name": atom_name,
                    "residue_name": line[17:20].strip(),
                    "chain_id": line[21].strip(),
                    "residue_seq": int(line[22:26].strip()),
                    "x": float(line[30:38].strip()),
                    "y": float(line[38:46].strip()),
                    "z": float(line[46:54].strip()),
                    "occupancy": float(line[54:60].strip()),
                    "temp_factor": float(line[60:66].strip()),
                    "element": line[76:78].strip()
                }
                atoms.append(atom_info)
            elif line.startswith("CONECT"):
                record = [int(x) for x in line.split()[1:]]
                conect_records.append(record)

    atoms_df = pd.DataFrame(atoms)
    atoms_data = atoms_df[['atom_name', 'element', 'x', 'y', 'z']].copy()
    atoms_data['index'] = range(1, len(atoms_data) + 1)

    neighbors_dict = {i+1: [] for i in range(len(atoms_df))}
    for record in conect_records:
        atom_id = record[0]
        bonded_atoms = record[1:]
        for bonded_atom in bonded_atoms:
            neighbors_dict[atom_id].append(bonded_atom)

    structured_data = []
    for i, row in atoms_data.iterrows():
        structured_data.append(f'''{row['index']}=["{row['atom_name']}", "{row['element']}", {row['x']}, {row['y']}, {row['z']}, {neighbors_dict[row['index']]}, {row['index']}]''')
    
    return structured_data

def load_mapping(file_path):
    mapping = {}
    with open(file_path, 'r') as file:
        for line in file:
            original, new = line.strip().split()
            mapping[original] = new
    return mapping

def find_link(x,structured_data):
    for element in structured_data:
        index = element[0]
        atom_name = element.split('"')[1]
        if atom_name == x:
            return index
    return None

def atom_count(structured_data):
    for element in structured_data:
        index = element.split('=')[0]
    return index

def print_header(clink, nlink, resolution, polymer, atomcount):

    print(f'res={atom_info['residue_name']}')
    print(f'resolution={resolution}')
    print(f'polymer={polymer}')
    print(f'atoms={atomcount}')
    print(f'link=[{clink}, {nlink}]')
    return

def main():
    parser = argparse.ArgumentParser(description='Process a PDB file with optional atom name mapping.')
    parser.add_argument('-f', '--file', type=str, required=True, help='Input PDB file')
    parser.add_argument('-m', '--mapping', type=str, help='Atom name mapping file')
    parser.add_argument('-n', '--nlink', type=str, help='Atom name for N link', default='N')
    parser.add_argument('-c', '--clink', type=str, help='Atom name for C link', default='C')
    parser.add_argument('-r', '--resolution', type=str, help='Model Resolution', default='atomistic')
    parser.add_argument('-p', '--polymer', type=int, help='Model Polymer', default=0)

    args = parser.parse_args()

    atom_mapping = {}
    if args.mapping:
        atom_mapping = load_mapping(args.mapping)
    
    structured_data = parse_pdb(args.file, atom_mapping)

    clink = find_link(args.clink, structured_data)
    nlink = find_link(args.nlink, structured_data)
    atomcount = atom_count(structured_data)

    print_header(clink, nlink, args.resolution, args.polymer, atomcount)
    
    for data in structured_data:
        print(data)

if __name__ == '__main__':
    main()
