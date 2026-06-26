import os
os.environ['MPLBACKEND'] = 'Agg'  # Hard-overrides the broken environment variable

# Now your original imports can follow safely below
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

def generate_decomposition_plot(csv_file, output_png, output_csv):
    if not os.path.exists(csv_file):
        print(f"❌ Error: {csv_file} not found.")
        return

    print(f"📂 Reading binding energy data from: {csv_file}...")
    
    # 1. Load the data 
    # Skipping the first few lines of metadata to reach the header
    try:
        df = pd.read_csv(csv_file, skiprows=3)
        # If 'Frame #' is not the first column, try skipping 4 lines instead
        if 'Frame #' not in df.columns:
            df = pd.read_csv(csv_file, skiprows=4)
    except Exception as e:
        print(f"❌ Error reading CSV: {e}")
        return

    # 2. Clean the data
    # Remove any non-numeric rows (like 'Receptor:' or text banners inside the file)
    df = df[pd.to_numeric(df['TOTAL'], errors='coerce').notnull()]
    
    # Convert energy columns to numeric
    energy_cols = ['van der Waals', 'Electrostatic', 'Polar Solvation', 'Non-Polar Solv.', 'TOTAL']
    for col in energy_cols:
        df[col] = pd.to_numeric(df[col])

    # 3. Calculate Statistics (Mean and Std Dev) across all frames
    summary = df.groupby('Residue')['TOTAL'].agg(['mean', 'std']).reset_index()
    
    # Sort by mean energy (from strongest favorable binders to least favorable)
    summary = summary.sort_values(by='mean')

    # 4. Plotting
    plt.figure(figsize=(20, 8))
    
    # Color coding: Green for favorable binding energy contributions, Red for unfavorable
    colors = ['#27ae60' if x < 0 else '#c0392b' for x in summary['mean']]
    
    # Create bar chart with error bars representing Standard Deviation
    plt.bar(summary['Residue'], summary['mean'], 
            yerr=summary['std'], 
            color=colors, 
            edgecolor='black', 
            alpha=0.8,
            error_kw={'capsize': 3, 'elinewidth': 0.6, 'ecolor': '#2c3e50'})

    # Chart Styling
    plt.axhline(0, color='black', linewidth=1.2)
    plt.ylabel('Binding Energy Contribution ($\Delta$G kcal/mol)', fontweight='bold', fontsize=12)
    plt.xlabel('Protein Residues (Chain:Residue:Number)', fontweight='bold', fontsize=12)
    
    # REVISION 1: Title customized to explicitly display the fxn protein-ligand footprint
    plt.title('Protein-Ligand (fxn) Per-Residue Energy Decomposition Summary', fontsize=16, pad=25, fontweight='bold')
    
    plt.xticks(rotation=90, fontsize=9)
    plt.grid(axis='y', linestyle='--', alpha=0.4)
    plt.tight_layout()

    # 5. Save Results directly back to Google Drive paths
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    summary.to_csv(output_csv, index=False)
    
    print("-" * 50)
    print(f"✅ Decomposition chart saved to: {output_png}")
    print(f"📊 Statistical summary saved to: {output_csv}")
    print("-" * 50)
    print("🌟 Top 5 Strongest Interacting Residues with Ligand (fxn):")
    print(summary.head(5)[['Residue', 'mean']])

if __name__ == "__main__":
    # REVISION 2: Hardcoded absolute Google Drive paths for your Protein-Ligand workspace
    base_dir = "/content/drive/MyDrive/ColabMD-Edu_Protein-Ligand/Analysis/MMGBSA/"
    
    target_csv = os.path.join(base_dir, "FINAL_DECOMP_MMPBSA.csv") 
    output_image = os.path.join(base_dir, "mmpbsa_decomposition_chart.png")
    output_summary_csv = os.path.join(base_dir, "residue_decomposition_summary.csv")
    
    generate_decomposition_plot(target_csv, output_image, output_summary_csv)