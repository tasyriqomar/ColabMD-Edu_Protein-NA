import os

def split_pdb(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    output_dir = os.path.join(os.path.dirname(file_path), 'split_pdbs')
    os.makedirs(output_dir, exist_ok=True)
    
    title = None
    pdb_lines = []
    pdb_index = 0

    for line in lines:
        if line.startswith("TITLE"):
            if pdb_lines:  # If there are accumulated lines, write them to a new file
                output_file = os.path.join(output_dir, f"pdb_{pdb_index}.pdb")
                with open(output_file, 'w') as out_file:
                    out_file.writelines([title] + pdb_lines)
                pdb_index += 1
                pdb_lines = []  # Reset the pdb_lines for the next PDB segment
            title = line  # Store the title line
        elif line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("TER"):
            pdb_lines.append(line)  # Accumulate PDB lines

    # Write the last PDB segment if exists
    if pdb_lines:
        output_file = os.path.join(output_dir, f"pdb_{pdb_index}.pdb")
        with open(output_file, 'w') as out_file:
            out_file.writelines([title] + pdb_lines)
            
    print(f"Successfully split the PDB file into {pdb_index + 1} individual files in the directory '{output_dir}'.")

if __name__ == "__main__":
    current_dir = os.getcwd()
    input_pdb_file = os.path.join(current_dir, "concatenated.pdb")  # Assuming the concatenated PDB file is in the current directory
    split_pdb(input_pdb_file)
