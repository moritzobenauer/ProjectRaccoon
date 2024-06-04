import numpy as np
from pynput import keyboard
import webbrowser
import os
import pyautogui
import argparse

def parse_pdb(file_path):
    atoms = []
    coords = []

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_info = line[12:16].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                atoms.append(atom_info)
                coords.append([x, y, z])

    return atoms, np.array(coords)

def coords_to_pdb(atoms, coords, original_lines):
    pdb_lines = []
    atom_index = 0
    for line in original_lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            coord = coords[atom_index]
            pdb_lines.append(f"{line[:30]}{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}{line[54:]}")
            atom_index += 1
        else:
            pdb_lines.append(line)
    return "".join(pdb_lines)

# Function to generate HTML content for visualization
def generate_html(pdb_strs):
    script = ""
    for pdb_str in pdb_strs:
        script += f"""
        viewer.addModel(`{pdb_str}`, "pdb");
        viewer.setStyle({{sphere: {{radius: 0.8}}}});
        """

    return f"""
<!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Project RACCOON: Geometry Builder</title>
        <script src="https://3dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
        <link href="https://fonts.googleapis.com/css2?family=Roboto:wght@400;700&display=swap" rel="stylesheet">
        <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css">
        <link rel="stylesheet" href="style.css">
    </head>
    
    <body>
        <div class=container>
        <header>
            <div id="header">
            <img src="logo.png" id="logo">
            <h2>Geometry Builder </h2>
            </div>           
            <p>Selected Molecule: {molecule_files[selected_molecule_index]}</p>
            
        </header>
        </div>
        <div id="viewer_wrapper">

        <div id="viewer"></div>
        <script>
            let viewer = $3Dmol.createViewer("viewer", {{backgroundColor: "white"}});
            {script}
            viewer.zoomTo();
            viewer.render();

            function refreshViewer() {{
                viewer.removeAllModels();
                {script}
                viewer.zoomTo();
                viewer.render();
            }}

            window.onbeforeunload = refreshViewer;
        </script>
        
        </div>
 

    <div class="container">
        
        <section class="content">
            <div style='text-align: left'>
                <h2>Translation & Rotation</h2>
                <p><kbd>W</kbd>, <kbd>A</kbd>, <kbd>S</kbd>, <kbd>D</kbd>: rotation</p>
                <p><kbd>Q</kbd>, <kbd>E</kbd>: mirroring</p>
                <p><kbd>U</kbd>, <kbd>N</kbd>: up and down</p>
                <p><kbd>H</kbd>, <kbd>J</kbd>: left and right</p>
                <p></p>
            </div>
            <div style='text-align: left'>
                <h2>File I/O</h2>
                <p><kbd>C</kbd>: Switch selected molecule</p>
                <p><kbd>1</kbd>: Exit without saving.</p>
                <p><kbd>2</kbd>: Exit with saving.</p>
            </div>
            <div style='text-align: justify'>
                <h2>About</h2>
                <p>This tool uses the <a href="https://github.com/avirshup/py3dmol" class="link">Py3Dmol</a> libary to display molecules. 
                Please reference Py3Dmol as well as <a href="https://github.com/moritzobenauer/ProjectRaccoon" class="link">Project RACCOON</a> when using this tool. Documentation available via github.</p>
            </div>
        </section>
        <footer>
            <i class="fa fa-github" style="font-size:15px"></i>
            <a href="https://github.com/moritzobenauer/ProjectRaccoon" class="link">projectraccoon.org </a>
            
            <i class="fa fa-book" style="font-size:15px; margin-left: 2px;"> </i>
            <a href="https://doi.org/10.21105/joss.06293" class="link">10.21105/joss.06293 </a>
             
        </footer>
    </div>
    </body>
    </html>
    """

def update_view():

    global atoms, coords, original_lines, selected_molecule_index
    atoms, coords = molecules[selected_molecule_index]
    original_lines = original_lines_list[selected_molecule_index]
    pdb_strs = [coords_to_pdb(molecules[i][0], molecules[i][1], original_lines_list[i]) for i in range(len(molecules))]
    html_content = generate_html(pdb_strs)
    with open("molecule_viewer.html", "w") as file:
        file.write(html_content)
    pyautogui.hotkey('f5')

def on_press(key):
    global coords, selected_molecule_index, save_results
    try:
        if key.char == 'w':
            rotation_matrix = np.array([[1, 0, 0],
                                        [0, np.cos(rotate_step), -np.sin(rotate_step)],
                                        [0, np.sin(rotate_step), np.cos(rotate_step)]])
            coords[:] = np.dot(coords, rotation_matrix)
        elif key.char == 's':
            rotation_matrix = np.array([[1, 0, 0],
                                        [0, np.cos(-rotate_step), -np.sin(-rotate_step)],
                                        [0, np.sin(-rotate_step), np.cos(-rotate_step)]])
            coords[:] = np.dot(coords, rotation_matrix)
        elif key.char == 'a':
            rotation_matrix = np.array([[np.cos(rotate_step), 0, np.sin(rotate_step)],
                                        [0, 1, 0],
                                        [-np.sin(rotate_step), 0, np.cos(rotate_step)]])
            coords[:] = np.dot(coords, rotation_matrix)
        elif key.char == 'd':
            rotation_matrix = np.array([[np.cos(-rotate_step), 0, np.sin(-rotate_step)],
                                        [0, 1, 0],
                                        [-np.sin(-rotate_step), 0, np.cos(-rotate_step)]])
            coords[:] = np.dot(coords, rotation_matrix)
        elif key.char == 'q':
            rotation_matrix = np.array([[np.cos(rotate_step), -np.sin(rotate_step), 0],
                                        [np.sin(rotate_step), np.cos(rotate_step), 0],
                                        [0, 0, 1]])
            coords[:] = np.dot(coords, rotation_matrix)
        elif key.char == 'e':
            rotation_matrix = np.array([[np.cos(-rotate_step), -np.sin(-rotate_step), 0],
                                        [np.sin(-rotate_step), np.cos(-rotate_step), 0],
                                        [0, 0, 1]])
            coords[:] = np.dot(coords, rotation_matrix)
        elif key.char == '2':
            print('Save & Exit')
            save_results = True
            return False
        elif key.char == '1':
            print('Exit')
            return False
        elif key.char == 'u':
            coords[:, 1] += move_step
        elif key.char == 'n':
            coords[:, 1] -= move_step
        elif key.char == 'h':
            coords[:, 0] -= move_step
        elif key.char == 'j':
            coords[:, 0] += move_step
        elif key.char == 'c':
            selected_molecule_index = (selected_molecule_index + 1) % len(molecules)
            print(selected_molecule_index)
        update_view()
    except AttributeError:
        pass

def main():

    global molecule_files, move_step, rotate_step, molecules, selected_molecule_index, original_lines, original_lines_list, sphere_size, save_results

    parser = argparse.ArgumentParser(description="")
    
    parser.add_argument('-f', '--files', nargs='+', help='PDB files', required=True)
    parser.add_argument('-m', '--trans', type=float, help='Translation Step Size', default=0.1)
    parser.add_argument('-r', '--rot', type=float, help='Rotation Step Size', default=5.0)
    parser.add_argument('-s', '--sphere', type=float, help='Sphere Size for py3Dmol', default=0.8)

    args = parser.parse_args()
    
    molecule_files = args.files
    move_step = args.trans
    rotate_step = np.radians(args.rot)
    

    molecules = [parse_pdb(file) for file in molecule_files]
    original_lines_list = [open(file).readlines() for file in molecule_files]
    selected_molecule_index = 0

    # Update atoms and coords to the first molecule initially
    atoms, coords = molecules[selected_molecule_index]
    original_lines = original_lines_list[selected_molecule_index]

    webbrowser.open("molecule_viewer.html", new=0)
    update_view()

    with keyboard.Listener(on_press=on_press) as listener:
        listener.join()

    if save_results == True:
        for i, (atoms, coords, original_lines) in enumerate(zip([mol[0] for mol in molecules], [mol[1] for mol in molecules], original_lines_list)):
            updated_pdb_str = coords_to_pdb(atoms, coords, original_lines)
            with open(molecule_files[i].replace('.pdb', '_updated.pdb'), 'w') as file:
                file.write(updated_pdb_str)

        print("Updated PDB files saved.")

    else:
        print('No PDB files saved.')
    
if __name__ == '__main__':
    save_results = False
    main()
