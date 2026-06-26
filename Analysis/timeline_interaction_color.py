import os
import pandas as pd
import matplotlib
# Force non-interactive backend to avoid ValueError in Conda/Colab
matplotlib.use('Agg') 
import matplotlib.pyplot as plt

# Function to safely parse frame integers from the pdbid string
def parse_frame_num(pdbid_str):
    try:
        # Splits 'PDB_89ns_fxn' -> '89ns' -> '89' -> 89
        return int(pdbid_str.split('_')[1].replace('ns', ''))
    except (IndexError, ValueError):
        return -1

# Function to extract frames from CSV files
def extract_frames(csv_files):
    frames = set()
    for csv_file in csv_files:
        df = pd.read_csv(csv_file)
        frames.update(df['pdbid'].apply(parse_frame_num))
    if -1 in frames:
        frames.remove(-1) # Clean out any unparseable rows
    return sorted(frames)

# Function to create timeline
def create_timeline(frames, csv_files):
    timeline = {}
    for frame in frames:
        timeline[frame] = {}
        for csv_file in csv_files:
            df = pd.read_csv(csv_file)
            
            # REVISION 1: Dynamic frame matching instead of hardcoded '_PROTEIN' string
            df['frame_idx'] = df['pdbid'].apply(parse_frame_num)
            
            if frame in df['frame_idx'].values:
                bond_type = os.path.splitext(os.path.basename(csv_file))[0]
                if bond_type not in timeline[frame]:
                    timeline[frame][bond_type] = 1
                else:
                    timeline[frame][bond_type] += 1
                    
        if not timeline[frame]:
            timeline[frame] = {'No Bond': 1}
    return timeline

# Main function
def main():
    # REVISION 2: Updated absolute path to point to your Protein-Ligand project directory
    base_dir = "/content/drive/MyDrive/ColabMD-Edu_Protein-NA/Analysis/split_pdbs/"
    bond_csv_folder = os.path.join(base_dir, "bond_csv_files")
    
    if not os.path.exists(bond_csv_folder):
        print(f"❌ Error: Folder '{bond_csv_folder}' not found. Please check your path.")
        return

    csv_files = [os.path.join(bond_csv_folder, file) for file in os.listdir(bond_csv_folder) if file.endswith('.csv')]

    if not csv_files:
        print(f"⚠️ No CSV files found inside: {bond_csv_folder}")
        return

    frames = extract_frames(csv_files)
    timeline = create_timeline(frames, csv_files)

    # Prepare data for stacked bar chart
    data = {bond_type: [] for bond_type in csv_files}
    for frame in timeline:
        for bond_type in data:
            data[bond_type].append(timeline[frame].get(os.path.splitext(os.path.basename(bond_type))[0], 0))

    # Plot stacked bar chart
    plt.figure(figsize=(15, 6))
    bottom = None
    colors = {
        'hydrophobic_interactions': 'grey',
        'hydrogen_bonds': 'royalblue',
        'salt_bridges': 'gold',
        'halogen_bonds': 'cyan',
        'pi_cation_interactions': 'maroon',
        'pi_stacks': 'forestgreen',
        'water_bridges': 'mediumpurple'
    }

    for bond_type in data:
        # Clean label for legend
        clean_label = bond_type.replace('_', ' ').title().split("/")[-1].split(".")[0]
        color_key = os.path.splitext(os.path.basename(bond_type))[0]
        
        if bottom is None:
            plt.bar(frames, data[bond_type], label=clean_label, color=colors.get(color_key, 'black'), alpha=0.8, edgecolor='black')
            bottom = data[bond_type]
        else:
            plt.bar(frames, data[bond_type], bottom=bottom, label=clean_label, color=colors.get(color_key, 'black'), alpha=0.8, edgecolor='black')
            bottom = [bottom[i] + data[bond_type][i] for i in range(len(bottom))]

    plt.xlabel('Simulation Time (ns)', fontweight='bold')
    plt.ylabel('Interaction Type Presence Count', fontweight='bold')
    
    # REVISION 3: Title customized for ligand analysis
    plt.title('Protein-Ligand (fxn) Non-Bonded Interaction Timeline', fontsize=14, pad=15, fontweight='bold')
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="Interaction Types")
    plt.grid(True, axis='y', linestyle='--', alpha=0.4)
    
    # REVISION 4: Dynamic X-axis ticks matching your actual frame array (e.g., 89 to 100)
    plt.xticks(frames)
    plt.tight_layout()

    # Save output plot back to your Analysis directory
    output_filename = os.path.join(base_dir, "bond_timeline_plot.png")
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"✅ Success: Timeline plot saved to: {output_filename}")
    plt.close()

if __name__ == "__main__":
    main()
