import os
import sys
import shutil
from Bio import PDB

# GROMACS-compatible residue mappings (AMBER / CHARMM force fields)
DNA_MAPPING = {
    "A": "DA", "T": "DT", "U": "DT", "G": "DG", "C": "DC",
    "DA": "DA", "DT": "DT", "DG": "DG", "DC": "DC",
    "RA": "DA", "RU": "DT", "RG": "DG", "RC": "DC",
    "ADE": "DA", "THY": "DT", "GUA": "DG", "CYT": "DC",
}

RNA_MAPPING = {
    "A": "A", "U": "U", "G": "G", "C": "C",
    "RA": "A", "RU": "U", "RG": "G", "RC": "C",
    "ADE": "A", "URA": "U", "GUA": "G", "CYT": "C",
}

ATOM_NAME_FIXES = {
    "O1P": "OP1",
    "O2P": "OP2",
}


def detect_nucleic_acid_type(structure) -> str:
    """
    Detects whether a structure is RNA or DNA based on Uracil residues
    or ribose 2'-hydroxyl (O2') atoms.
    """
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.resname.strip().upper()
                if resname in {"U", "RU", "URA"}:
                    return "RNA"
                for atom in residue:
                    atom_name = atom.name.strip().replace("*", "'").upper()
                    if atom_name in {"O2'", "HO2'"}:
                        return "RNA"
    return "DNA"


def standardize_atom_name(atom):
    """Converts legacy asterisk (*) sugar notation and O1P/O2P to GROMACS/PDBv3 names."""
    clean_name = atom.name.strip().replace("*", "'")
    clean_name = ATOM_NAME_FIXES.get(clean_name, clean_name)
    atom.name = clean_name
    atom.fullname = f" {clean_name:<3}" if len(clean_name) < 4 else f"{clean_name:<4}"


# 1. Determine input file (supports both aptamer_RNA.pdb and aptamer_DNA.pdb)
if os.path.exists("aptamer_RNA.pdb"):
    input_pdb = "aptamer_RNA.pdb"
    explicit_type = "RNA"
elif os.path.exists("aptamer_DNA.pdb"):
    input_pdb = "aptamer_DNA.pdb"
    explicit_type = None
else:
    raise FileNotFoundError("Neither aptamer_RNA.pdb nor aptamer_DNA.pdb was found.")

# Allow optional CLI override: `python rename_dna.py RNA` or `python rename_dna.py DNA`
if len(sys.argv) > 1 and sys.argv[1].upper() in {"DNA", "RNA"}:
    explicit_type = sys.argv[1].upper()

parser = PDB.PDBParser(QUIET=True)
structure = parser.get_structure("apt", input_pdb)

# 2. Resolve nucleic acid type (explicit flag or auto-detect from structure)
aptamer_type = explicit_type or detect_nucleic_acid_type(structure)

if aptamer_type == "RNA" and not os.path.exists("aptamer_RNA.pdb"):
    shutil.copy(input_pdb, "aptamer_RNA.pdb")

mapping = RNA_MAPPING if aptamer_type == "RNA" else DNA_MAPPING

# 3. Standardize residue names and sugar/phosphate atom names
for model in structure:
    for chain in model:
        for residue in chain:
            resname = residue.resname.strip().upper()
            if resname in mapping:
                residue.resname = mapping[resname]

            for atom in list(residue):
                standardize_atom_name(atom)
                # Strip 2'-hydroxyl atoms only when targeting DNA
                if aptamer_type == "DNA" and atom.name.strip() in {"O2'", "HO2'"}:
                    residue.detach_child(atom.id)

io = PDB.PDBIO()
io.set_structure(structure)
io.save("aptamer_fixed.pdb")

print(f"Saved fixed {aptamer_type} PDB: aptamer_fixed.pdb")
