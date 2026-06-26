# 1. Navigate to your directory
cd "/content/drive/MyDrive/ColabMD-Edu_Protein-Ligand/Analysis/split_pdbs/"

# 2. Loop through every PDB file
for pdb in pdb_*.pdb; do
    # Get the name without the .pdb extension (e.g., pdb_0)
    base_name="${pdb%.pdb}"
    
    # Create the folder if it doesn't exist
    if [ ! -d "$base_name" ]; then
        mkdir -p "$base_name"
    fi
    
    echo "------------------------------------------"
    echo "Analyzing $pdb --> Folder: $base_name"
    
    # Run PLIP
    # -f: input file
    # -o: output directory (now matches the pdb name)
    # -y: generate PyMOL session
    # -x: generate XML report
 
    plip -f "$pdb" -p -y -x -t -o "$base_name"
done

echo "------------------------------------------"
echo "Analysis Complete. Folders created for pdb_0 through pdb_10."