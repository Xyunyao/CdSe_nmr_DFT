import os
import sys
import nmrglue as ng
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Add the path to the directory containing scripts_module
sys.path.append(os.path.abspath('/Users/yunyao_1/Dropbox/CdSe-data-owen/scripts_module'))

from orca_analyzer import ShieldingTensorAnalyzer
from read_orca import extract_tensors_from_file

def main(path):
    # Extract tensors for Cd and Se
    output_dir = os.path.dirname(path)
    file_name = os.path.splitext(os.path.basename(path))[0] 
    cd_dataset = extract_tensors_from_file(path, 'Cd')
    se_dataset = extract_tensors_from_file(path, 'Se')
    
    # Initialize analyzers
    cd_analyzer = ShieldingTensorAnalyzer(cd_dataset, 3059.972)
    se_analyzer = ShieldingTensorAnalyzer(se_dataset, 836.688)

    # Compute isotropic shift, delta value, and ita value for Cd
    cd_iso = [cd_analyzer.isotropic_shift(key) for key in cd_dataset.keys() if 'Cd' in key]
    cd_delta = [cd_analyzer.delta_value(key) for key in cd_dataset.keys() if 'Cd' in key]
    cd_ita = [cd_analyzer.ita_value(key) for key in cd_dataset.keys() if 'Cd' in key]

    # Compute isotropic shift, delta value, and ita value for Se
    se_iso = [se_analyzer.isotropic_shift(key) for key in se_dataset.keys() if 'Se' in key]
    se_delta = [se_analyzer.delta_value(key) for key in se_dataset.keys() if 'Se' in key]
    se_ita = [se_analyzer.ita_value(key) for key in se_dataset.keys() if 'Se' in key]

    # Save cd_iso values and their indices to a CSV file
    cd_indices = list(range(len(cd_iso)))  # Generate indices for Cd isotropic shifts
    cd_iso_df = pd.DataFrame({
        'Cd Index': cd_indices,
        'Cd Isotropic Shift': cd_iso
    })
    csv_path = os.path.join(output_dir, f'{file_name}_cd_iso_data.csv')
    cd_iso_df.to_csv(csv_path, index=False)  # Save the dataframe to CSV

    # Plot Cd data
    x = np.arange(len(cd_iso))
    fig_cd, axs_cd = plt.subplots(3, 1, figsize=(8, 12))

    axs_cd[0].bar(x, cd_iso, color='blue')
    axs_cd[0].set_title('cd_iso')
    axs_cd[0].set_xlabel('Index')
    axs_cd[0].set_ylabel('Values')

    axs_cd[1].bar(x, cd_ita, color='orange')
    axs_cd[1].set_title('cd_ita')
    axs_cd[1].set_xlabel('Index')
    axs_cd[1].set_ylabel('Values')

    axs_cd[2].bar(x, cd_delta, color='green')
    axs_cd[2].set_title('cd_delta')
    axs_cd[2].set_xlabel('Index')
    axs_cd[2].set_ylabel('Values')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{file_name}_Cd_plot.png'))  # Save Cd plot
    plt.show()

    # Plot Se data
    x = np.arange(len(se_iso))
    fig_se, axs_se = plt.subplots(3, 1, figsize=(8, 12))

    axs_se[0].bar(x, se_iso, color='blue')
    axs_se[0].set_title('se_iso')
    axs_se[0].set_xlabel('Index')
    axs_se[0].set_ylabel('Values')

    axs_se[1].bar(x, se_ita, color='orange')
    axs_se[1].set_title('se_ita')
    axs_se[1].set_xlabel('Index')
    axs_se[1].set_ylabel('Values')

    axs_se[2].bar(x, se_delta, color='green')
    axs_se[2].set_title('se_delta')
    axs_se[2].set_xlabel('Index')
    axs_se[2].set_ylabel('Values')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{file_name}_Se_plot.png'))  # Save Se plot
    plt.show()

    # Function to generate a Gaussian distribution
    def gaussian(x, scale, mu, sigma):
        return scale * np.exp(-0.5 * ((x - mu) / sigma)**2)

    # Example data
    group1 = [cd_iso[3], cd_iso[6], cd_iso[7], cd_iso[11]]
    group2 = [x for x in cd_iso if x not in group1]

    # Create a figure and axis
    fig=plt.figure(figsize=(8, 6))
    ax = plt.gca()

    # Plot the histogram and Gaussian fit for each group
    hist1, bins1 = np.histogram(group1, bins=20, density=False)
    plt.hist(group1, bins=10, density=False, alpha=0.5, color='blue', label='Face Cd')

    hist2, bins2 = np.histogram(group2, bins=20, density=False)
    plt.hist(group2, bins=10, density=False, alpha=0.5, color='orange', label='Edge Cd')

    # Parameters for the Gaussian distributions
    params1 = [max(hist1), np.mean(group1), np.std(group1)]
    params2 = [max(hist2), np.mean(group2), np.std(group2)]

    # Plot the Gaussian distribution for each group
    x = np.linspace(min(bins1[:-1]), max(bins1[:-1]), 100)
    plt.plot(x, gaussian(x, *params1), 'b-',label=f'Face Cd ($\mu={params1[1]:.2f}$, $\sigma={params1[2]:.2f}$)')

    x = np.linspace(min(bins2[:-1]), max(bins2[:-1]), 100)
    plt.plot(x, gaussian(x, *params2), 'r-',label=f'Edge Cd ($\mu={params2[1]:.2f}$, $\sigma={params2[2]:.2f}$)')

    # Set labels and title
    plt.xlabel('Calculated $^{113}$Cd chemical shift')
    plt.ylabel('Frequency')
    plt.title('Gaussian Distributions for Each Group')

    # Add legend
    plt.legend()
    plt.savefig(os.path.join(output_dir, f'{file_name}_Cd_distribution_plot.png'))
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python QM_process.py <path>")
        sys.exit(1)

    path = sys.argv[1]
    main(path)

