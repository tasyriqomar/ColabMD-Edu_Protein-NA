import json
import csv
import os
import glob
import matplotlib
# FIX: Use 'Agg' for Colab/Conda environments
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from collections import defaultdict, Counter

def extract_residues(json_data):
    """ Extracts interacting residues and organizes them by general type and specific chain. """
    interaction_residues = defaultdict(list)
    chain_residues = defaultdict(lambda: defaultdict(list))

    for pdb_entry in json_data.values():
        for binding_site in pdb_entry.get("bindingsites", []):
            interactions = binding_site.get("interactions", {})

            for interaction_type, bonds in interactions.items():
                for bond in bonds:
                    # Extracts chain, residue number, and type
                    chain = bond.get("reschain", "Unknown")
                    resnr = bond.get("resnr", "0")
                    restype = bond.get("restype", "UNK")
                    
                    formatted_residue = f"{restype}:{resnr}" # Format for the X-axis
                    
                    # Store for overall stats
                    interaction_residues[interaction_type].append(f"{chain}:{formatted_residue}")
                    
                    # Store specifically for that chain's plot
                    chain_residues[chain][interaction_type].append(formatted_residue)

    return interaction_residues, chain_residues

def save_chain_plot(chain_data, chain_id, output_dir):
    """ Generates and saves a bar plot for a specific chain. """
    if not chain_data:
        return

    color_map = {
        'hydrophobic_interactions': 'grey',
        'hydrogen_bonds': 'royalblue',
        'salt_bridges': 'orange',
        'pi_stacks': 'forestgreen',
        'water_bridges': 'mediumpurple',
        'pi_cation_interactions': 'maroon'
    }

    plt.figure(figsize=(12, 6))
    
    # Track labels for legend
    seen_labels = set()
    has_data = False

    for interaction_type, residues in chain_data.items():
        # Get top 8 residues for this specific interaction type on this chain
        top_res = Counter(residues).most_common(8)
        if not top_res:
            continue
            
        has_data = True
        color = color_map.get(interaction_type, 'black')
        
        for residue, count in top_res:
            label = interaction_type if interaction_type not in seen_labels else ""
            plt.bar(residue, count, color=color, edgecolor='black', label=label, alpha=0.8)
            seen_labels.add(interaction_type)

    if not has_data:
        plt.close()
        return

    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="Interaction Types")
    plt.xlabel(f'Protein Residues (Chain {chain_id})', fontweight='bold')
    plt.ylabel('Interaction Frequency (Across Frames)', fontweight='bold')
    
    # REVISION 1: Title customized to display ligand (fxn) context
    plt.title(f'Top Interacting Residues: Chain {chain_id} ↔ Ligand (fxn)', fontsize=14, pad=15)
    
    plt.xticks(rotation=45, ha="right")
    plt.grid(axis='y', linestyle='--', alpha=0.3)
    plt.tight_layout()
    
    output_path = os.path.join(output_dir, f"top_residues_chain_{chain_id}.png")
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"✅ Generated plot for Chain {chain_id}: {output_path}")

if __name__ == "__main__":
    # REVISION 2: Path updated to point to the Protein-Ligand project directory
    base_dir = "/content/drive/MyDrive/ColabMD-Edu_Protein-NA/Analysis/split_pdbs/"
    
    json_file = os.path.join(base_dir, "extracted_data.json")
    output_dir = os.path.join(base_dir, 'chain_analysis_results')
    os.makedirs(output_dir, exist_ok=True)

    if not os.path.exists(json_file):
        print(f"❌ Error: {json_file} not found!")
    else:
        with open(json_file, "r") as file:
            data = json.load(file)

        # 1. Extract data for all chains
        _, chain_data_map = extract_residues(data)

        # 2. Write a detailed CSV for all chains
        csv_path = os.path.join(output_dir, "detailed_residue_counts_by_chain.csv")
        with open(csv_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Chain', 'Interaction_Type', 'Residue', 'Count'])
            
            for chain, interactions in sorted(chain_data_map.items()):
                for itype, res_list in interactions.items():
                    counts = Counter(res_list)
                    for res, count in counts.most_common():
                        writer.writerow([chain, itype, res, count])
        
        print(f"📊 CSV report generated: {csv_path}")

        # 3. Generate a plot for EVERY chain detected mapping to fxn
        print("\n🎨 Generating visualizations...")
        for chain_id in sorted(chain_data_map.keys()):
            save_chain_plot(chain_data_map[chain_id], chain_id, output_dir)

        print(f"\n🚀 Analysis complete. All results are in: {output_dir}")
