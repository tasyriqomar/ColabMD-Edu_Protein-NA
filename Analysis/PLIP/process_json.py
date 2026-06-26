import json
import csv
import os

# Read JSON data from file
with open('extracted_data.json', 'r') as file:
    data = json.load(file)

# Initialize a dictionary to hold bond type data
bond_data = {}

# Iterate over each PDB entry
for pdb_entry in data.values():
    pdb_id = pdb_entry["pdbid"]

    # Iterate through all binding sites
    for binding_site in pdb_entry["bindingsites"]:
        interactions = binding_site["interactions"]

        for bond_type, bonds in interactions.items():
            if bonds:  # Check if there are any interactions of this type
                if bond_type not in bond_data:
                    bond_data[bond_type] = []

                for bond in bonds:
                    bond["pdbid"] = pdb_id  # Add the PDB ID to each bond entry
                    bond["chain"] = binding_site["identifiers"]["chain"]  # Include chain info
                    bond_data[bond_type].append(bond)

# Create a directory to store the CSV files
output_dir = 'bond_csv_files'
os.makedirs(output_dir, exist_ok=True)

# Generate CSV files for each bond type
for bond_type, bonds in bond_data.items():
    if bonds:  # Check if there are any bonds to write
        csv_file = os.path.join(output_dir, f"{bond_type}.csv")
        with open(csv_file, mode='w', newline='') as file:
            fieldnames = bonds[0].keys()
            writer = csv.DictWriter(file, fieldnames=fieldnames)
            writer.writeheader()
            for bond in bonds:
                writer.writerow(bond)
        print(f"Generated CSV file for {bond_type}: {csv_file}")
