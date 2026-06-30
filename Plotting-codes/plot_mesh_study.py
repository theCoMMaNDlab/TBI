# Author: Fatemeh Jafarabadi, Apr 2026
# Units: Length [micrometers], Time [days], Stress [Pa]
# ============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import re

# -----------------------------------------------------
# Configuration
# -----------------------------------------------------
TARGET_TIME = 0.25
RECOMMENDED_ELEMENTS = 18781
NT11_CSV = 'NT11_all.csv'
NT12_CSV = 'NT12_all.csv'
UMAG_CSV = 'UMAG_all.csv'
OUTPUT_PNG = 'mesh-t=025.png'

# -----------------------------------------------------
# Data Loading and Preparation
# -----------------------------------------------------
def load_and_prepare_data():
    nt11 = pd.read_csv(NT11_CSV)
    nt12 = pd.read_csv(NT12_CSV)
    umag = pd.read_csv(UMAG_CSV)

    time_idx = np.argmin(np.abs(nt11['Time'].values - TARGET_TIME))

    meshes = [col for col in nt11.columns if col != 'Time']
    element_counts = [int(re.search(r'elm=(\d+)', mesh_name).group(1)) for mesh_name in meshes]
    order = np.argsort(element_counts)
    element_counts = [element_counts[i] for i in order]

    pressure = nt11.iloc[time_idx, 1:].astype(float).values
    astrocytes = nt12.iloc[time_idx, 1:].astype(float).values
    displacement = umag.iloc[time_idx, 1:].astype(float).values

    pressure = pressure[order]
    astrocytes = astrocytes[order]
    displacement = displacement[order]

    pressure = pressure/float(pressure[0])
    astrocytes = astrocytes/float(astrocytes[0])
    displacement = displacement/float(displacement[0])

    return element_counts, displacement, pressure, astrocytes

# -----------------------------------------------------
# Plotting
# -----------------------------------------------------
def plot_output_vs_elements(x_vals, u_vals, p_vals, a_vals, output_name):
    plt.figure(figsize=(8, 6))
    plt.plot(x_vals, u_vals, '-o', color='black', label=r'$\tilde{u}$', linewidth=1.8)
    plt.plot(x_vals, p_vals, '-o', color='tab:blue', label=r'$\tilde{P}$', linewidth=1.8)
    plt.plot(x_vals, a_vals, '-o', color='tab:red', label=r'$\tilde{a}$', linewidth=1.8)

    plt.xlabel('Number of elements', fontsize=14)
    plt.ylabel('Normalized value', fontsize=14)
    plt.ylim(0.45, 1.55)
    plt.xlim(-500, 56500)
    plt.yticks([0.5, 0.75, 1.0, 1.25, 1.5])

    ax = plt.gca()
    ax.axvline(RECOMMENDED_ELEMENTS, color='green', linestyle='--',
               linewidth=1.2, label=f'{RECOMMENDED_ELEMENTS} elements')

    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax.tick_params(axis='y', which='minor', length=2, width=0.8)
    plt.legend(fontsize=14)
    plt.tick_params(axis='both', labelsize=14)
    plt.tight_layout()
    plt.savefig(output_name, dpi=200)

# -----------------------------------------------------
# Main Execution
# -----------------------------------------------------
def main():
    element_counts, displacement, pressure, astrocytes = load_and_prepare_data()
    plot_output_vs_elements(element_counts, displacement, pressure, astrocytes, OUTPUT_PNG)
    plt.show()


if __name__ == '__main__':
    main()
