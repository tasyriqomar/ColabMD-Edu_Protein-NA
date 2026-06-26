import os
os.environ['MPLBACKEND'] = 'Agg'  # Forces the headless file-saver backend

import numpy as np
import matplotlib
# CRITICAL: Use the 'Agg' backend to avoid display errors and enable file saving
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Function to read GROMACS .xvg files
def read_xvg(filename):
    data = []
    if not os.path.exists(filename):
        print(f"⚠️ Warning: {filename} not found.")
        return None
    with open(filename, 'r') as file:
        for line in file:
            # Skip comments (#) and metadata (@, &)
            if line.strip() and not line.startswith(('#', '@', '&')):
                try:
                    data.append(list(map(float, line.split())))
                except ValueError:
                    continue
    return np.array(data)

def main():
    # File names
    proj1_2_file = '2dproj1-2.xvg'
    proj1_3_file = '2dproj1-3.xvg'
    proj2_3_file = '2dproj2-3.xvg'

    # Read data
    proj1_2 = read_xvg(proj1_2_file)
    proj1_3 = read_xvg(proj1_3_file)
    proj2_3 = read_xvg(proj2_3_file)

    # Validate data presence
    if proj1_2 is None or proj1_3 is None or proj2_3 is None:
        print("❌ Error: Missing .xvg files. Please check your folder.")
        return

    # Plot Configuration (1000x800 equivalent: 10x8 inches at 100 DPI)
    plt.figure(figsize=(10, 8))

    # Scatter plots with hollow markers (as in your original script)
    plt.scatter(proj1_2[:, 0], proj1_2[:, 1], color='none', edgecolor='blue',
                marker='^', label='PC1 vs PC2', alpha=0.6)
    plt.scatter(proj1_3[:, 0], proj1_3[:, 1], color='none', edgecolor='red',
                marker='s', label='PC1 vs PC3', alpha=0.6)
    plt.scatter(proj2_3[:, 0], proj2_3[:, 1], color='none', edgecolor='green',
                marker='o', label='PC2 vs PC3', alpha=0.6)

    # Formatting
    plt.xlabel('Principal Component Projection (nm)', fontweight='bold')
    plt.ylabel('Principal Component Projection (nm)', fontweight='bold')
    plt.title('PCA Essential Dynamics: Combined 2D Projections', fontsize=14, pad=15)
    plt.legend(loc='best', frameon=True, shadow=True)
    plt.grid(True, linestyle='--', alpha=0.4)

    plt.tight_layout()

    # SAVE TO CURRENT DIRECTORY
    output_filename = "pca_projections_combined.png"
    plt.savefig(output_filename, dpi=300)
    plt.close()

    print(f"🚀 Success! Plot saved as: {os.getcwd()}/{output_filename}")

if __name__ == "__main__":
    main()