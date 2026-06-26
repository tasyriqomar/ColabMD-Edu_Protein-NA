import json
import os
import matplotlib
# Use 'Agg' backend for Colab/Conda compatibility
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from collections import defaultdict

def extract_interaction_types(json_data):
    """
    Extracts interaction counts from the structured JSON data.
    """
    interaction_types = defaultdict(int)

    for pdb_entry in json_data.values():
        for binding_site in pdb_entry.get("bindingsites", []):
            interactions = binding_site.get("interactions", {})
            for interaction_type, bonds in interactions.items():
                # Bonds is a list of interaction dictionaries
                interaction_types[interaction_type] += len(bonds)

    return interaction_types

def plot_interactions(interaction_counts, output_path):
    """
    Plots a bar chart and saves it to a file.
    """
    if not interaction_counts:
        print("⚠️ No interactions found in JSON. Check your extraction step.")
        return

    # Clean up names for the X-axis (e.g., remove underscores)
    types = list(interaction_counts.keys())
    clean_types = [t.replace('_', ' ').title() for t in types]
    counts = [interaction_counts[t] for t in types]

    # Color mapping
    colors_map = {
        'hydrophobic_interactions': 'grey',
        'hydrogen_bonds': 'royalblue',
        'salt_bridges': 'gold',
        'halogen_bonds': 'cyan',
        'pi_stacks': 'forestgreen',
        'pi_cation_interactions': 'maroon',
        'water_bridges': 'mediumpurple',
        'metal_complexes': 'crimson'
    }

    bar_colors = [colors_map.get(t, 'black') for t in types]

    plt.figure(figsize=(12, 7))
    bars = plt.bar(clean_types, counts, color=bar_colors, edgecolor='black', alpha=0.8)

    # Add count labels on top of bars
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, yval + 0.1, int(yval), ha='center', va='bottom', fontweight='bold')

    plt.xlabel('Interaction Types', fontweight='bold')
    plt.ylabel('Total Occurrences (All Frames)', fontweight='bold')
    
    # REVISION 1: Updated Title for Protein-Ligand (fxn)
    plt.title('Protein-Ligand (fxn) Interaction Profile Summary', fontsize=14, pad=20)
    
    plt.xticks(rotation=30, ha="right")
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.tight_layout()
    
    # Save the figure
    plt.savefig(output_path, dpi=300)
    print(f"✅ Bar chart saved to: {output_path}")
    plt.close()

if __name__ == "__main__":
    # REVISION 2: Updated Google Drive directory path to Protein-Ligand project folder
    base_dir = "/content/drive/MyDrive/ColabMD-Edu_Protein-Ligand/Analysis/split_pdbs/"
    
    json_file = os.path.join(base_dir, "extracted_data.json")
    output_img = os.path.join(base_dir, "interaction_summary_plot.png")

    if not os.path.exists(json_file):
        print(f"❌ Error: {json_file} not found!")
    else:
        with open(json_file, "r") as file:
            data = json.load(file)

        counts = extract_interaction_types(data)
        plot_interactions(counts, output_img)