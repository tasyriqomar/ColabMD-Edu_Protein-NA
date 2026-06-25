from Bio import PDB

mapping = {"A": "DA", "T": "DT", "G": "DG", "C": "DC"}

parser = PDB.PDBParser(QUIET=True)
structure = parser.get_structure("apt", "aptamer_DNA.pdb")

for model in structure:
    for chain in model:
        for residue in chain:
            resname = residue.resname.strip()
            if resname in mapping:
                residue.resname = mapping[resname]

io = PDB.PDBIO()
io.set_structure(structure)
io.save("aptamer_fixed.pdb")

print("Saved fixed PDB: aptamer_fixed.pdb")

