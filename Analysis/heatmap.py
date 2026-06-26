import os
# FIX 1: Explicitly force the non-interactive backend BEFORE importing matplotlib
os.environ['MPLBACKEND'] = 'Agg'

import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def get_section_data(lines, csv_file, start_marker, num_frames=11, num_residues=541):
    """Finds a section in the text and reads the numeric table via pandas."""
    start_idx = -1
    for i, line in enumerate(lines):
        if start_marker in line:
            # Find the header row specifically for this section
            for j in range(i, len(lines)):
                if "Frame #" in lines[j] and "Residue" in lines[j]:
                    start_idx = j
                    break
            if start_idx != -1: 
                break
    
    if start_idx == -1:
        return None
    
    # Read the data rows
    nrows = num_frames * num_residues
    df = pd.read_csv(csv_file, skiprows=start_idx, nrows=nrows)
    
    # Clean column names and use position to be safe
    df.columns = df.columns.str.replace('"', '').str.strip()
    
    # We strictly need Col 0 (Frame), Col 1 (Residue), and Col 7 (TOTAL)
    df = df.iloc[:, [0, 1, 7]]
    df.columns = ['Frame', 'Residue', 'Energy']
    
    # Convert to numeric
    df['Energy'] = pd.to_numeric(df['Energy'], errors='coerce')
    df['Frame'] = pd.to_numeric(df['Frame'], errors='coerce', downcast='integer')
    return df.dropna()

def generate_decomposition_heatmap(csv_file, output_png, num_frames=11, num_residues=541):
    if not os.path.exists(csv_file):
        print(f"❌ Error: {csv_file} not found.")
        return

    with open(csv_file, 'r') as f:
        lines = f.readlines()

    print(f"🔍 Parsing Complex, Receptor, and Ligand sections ({num_frames} frames, {num_residues} residues)...")
    comp = get_section_data(lines, csv_file, "Complex:", num_frames, num_residues)
    rec  = get_section_data(lines, csv_file, "Receptor:", num_frames, num_residues)
    lig  = get_section_data(lines, csv_file, "Ligand:", num_frames, num_residues)

    if comp is None or rec is None or lig is None:
        print("❌ Error: Could not isolate one of the sections. Check marker names inside the CSV.")
        return

    # 1. Calculate Binding Energy (Delta G)
    # G_bind = G_complex - (G_receptor + G_ligand)
    final_df = comp.copy()
    final_df['Binding_G'] = comp['Energy'] - (rec['Energy'] + lig['Energy'])

    # 2. Identify Top 20 Hotspots (lowest average energy / strongest binding interaction)
    top_20 = final_df.groupby('Residue')['Binding_G'].mean().sort_values().head(20).index
    plot_df = final_df[final_df['Residue'].isin(top_20)]

    # 3. Create Pivot Table for Heatmap
    pivot = plot_df.pivot_table(index='Residue', columns='Frame', values='Binding_G', aggfunc=np.mean)
    pivot = pivot.reindex(top_20)

    # 4. Visualization
    plt.figure(figsize=(14, 10))
    sns.heatmap(pivot, annot=True, fmt=".2f", cmap="RdYlGn_r", 
                linewidths=.5, cbar_kws={'label': 'Delta G Binding Energy Contribution (kcal/mol)'})
    
    # REVISION 1: Updated title and labels for the fxn protein-ligand system
    plt.title('MMPBSA Binding Energy Heatmap: Top 20 Hotspots Across Frames\nProtein-Ligand (fxn) System', fontsize=15, pad=20, fontweight='bold')
    plt.xlabel('Frame Number', fontweight='bold')
    plt.ylabel('Protein Residues (Chain:Type:Number)', fontweight='bold')
    
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    print(f"🚀 Success! Heatmap saved as: {output_png}")
    plt.close()

if __name__ == "__main__":
    # REVISION 2: Configured path routing for your absolute Google Drive workspace
    base_dir = "/content/drive/MyDrive/ColabMD-Edu_Protein-NA/Analysis/MMGBSA/"
    
    target_csv = os.path.join(base_dir, "FINAL_DECOMP_MMPBSA.csv")
    output_img = os.path.join(base_dir, "mmpbsa_heatmap_results.png")
    
    # REVISION 3: Defaulting to 541 residues to map your entire receptor structure automatically
    generate_decomposition_heatmap(target_csv, output_img, num_frames=11, num_residues=541)
