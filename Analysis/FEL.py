import os
os.environ['MPLBACKEND'] = 'Agg'  # Forces the headless file-saver backend

import numpy as np
import matplotlib
# Use 'Agg' for Colab/Conda compatibility
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_fel(filename, title, xlabel='PC2', ylabel='PC1', zlabel='G (kJ/mol)'):
    if not os.path.exists(filename):
        print(f"❌ Error: {filename} not found.")
        return

    # Load data
    # GROMACS sham output usually has 3 columns: PC1, PC2, Energy
    data = np.loadtxt(filename)

    x = data[:, 0]  # PC2
    y = data[:, 1]  # PC1
    z = data[:, 2]  # Gibbs energy

    # Create plot
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection='3d')

    # Surface plot with 'rainbow' or 'jet' (common for FEL)
    # antialiased=True makes the surface look smoother
    surf = ax.plot_trisurf(x, y, z, cmap='jet', edgecolor='none', alpha=0.9, antialiased=True)

    # Labeling
    ax.set_xlabel(xlabel, fontweight='bold')
    ax.set_ylabel(ylabel, fontweight='bold')
    ax.set_zlabel(zlabel, fontweight='bold')

    # Adjust view angle for better perspective of the energy wells
    ax.view_init(elev=30, azim=45)

    # Remove grid for a cleaner look
    ax.grid(False)

    # Color bar
    cbar = fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, pad=0.1)
    cbar.set_label(zlabel, fontweight='bold')

    plt.title(f'Free Energy Landscape: {title}', fontsize=15, pad=20)

    # SAVE TO CURRENT DIRECTORY
    output_png = "FEL_3D_plot.png"
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"🚀 Success! FEL plot saved as: {os.getcwd()}/{output_png}")

if __name__ == "__main__":
    # Generate plot for FEL.dat
    plot_fel('FEL.dat', 'SERT-fxn')