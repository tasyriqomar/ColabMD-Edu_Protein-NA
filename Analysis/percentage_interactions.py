import json
import os
import matplotlib
# Force non-interactive backend for stable cloud plotting
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from collections import defaultdict

def calculate_json_percentages(json_data):
    """
    Parses the structured JSON data, counts interaction frequencies 
    by counting bond list lengths, and calculates relative percentages.
    """
    total_counts = defaultdict(int)
    overall_total = 0

    # Dig into the JSON tree architecture
    for pdb_entry in json_data.values():
        for binding_site in pdb_entry.get("bindingsites", []):
            interactions = binding_site.get("interactions", {})
            for interaction_type, bonds in interactions.items():
                # Each 'bonds' entry is a list of interaction dictionaries
                bond_count = len(bonds)
                total_counts[interaction_type] += bond_count
                overall_total += bond_count

    if overall_total == 0:
        return {}

    # Calculate percentages based on total bond frequencies
    percentages = {
        interaction_type: (count / overall_total) * 100 
        for interaction_type, count in total_counts.items()
    }
    
    return percentages

def save_percentage_table(percentages, output_path):
    """ Renders a clean summary table image and saves it to Google Drive. """
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.axis('tight')
    ax.axis('off')
    
    # Title optimized for your fxn protein-ligand system
    plt.title('Protein-Ligand (fxn) Interaction Share Summary', fontsize=14, pad=10, fontweight='bold')
    
    table_data = [["Interaction Type", "Percentage Share"]]
    
    # Sort interactions dynamically from highest percentage to lowest
    for interaction_type, percentage in sorted(percentages.items(), key=lambda x: x[1], reverse=True):
        clean_name = interaction_type.replace('_', ' ').title()
        table_data.append([clean_name, f"{percentage:.2f}%"])
    
    table = ax.table(cellText=table_data, colLabels=None, cellLoc='center', loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1.2, 1.4)
    
    # Apply clean bold headers and subtle shading to the top row
    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_text_props(weight='bold')
            cell.set_facecolor('#f2f2f2')
            
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"✅ Summary table figure successfully saved to: {output_path}")

if __name__ == "__main__":
    # Absolute paths pointing directly to your active Protein-Ligand workspace
    base_dir = "/content/drive/MyDrive/ColabMD-Edu_Protein-NA/Analysis/split_pdbs/"
    json_file = os.path.join(base_dir, "extracted_data.json")
    output_img = os.path.join(base_dir, "interaction_percentages.png")
    
    if not os.path.exists(json_file):
        print(f"❌ Error: Target file '{json_file}' does not exist.")
    else:
        with open(json_file, 'r') as file:
            data = json.load(file)
            
        percentages = calculate_json_percentages(data)
        
        if not percentages:
            print("⚠️ Warning: No interactions were found inside your extracted_data.json.")
        else:
            save_percentage_table(percentages, output_img)
