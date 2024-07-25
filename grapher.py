import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import json
import tkinter as tk

from tkinter import filedialog


# Get filepath of json
# Create a root window and hide it
root = tk.Tk()
root.withdraw()

# Get the current working directory
current_directory = os.getcwd()

# Open file dialog with the current directory as the default
filepath = filedialog.askopenfilename(initialdir=current_directory, filetypes=[("JSON files", "*.json")])

# Load dictionary into 'data_file' variable
with open(filepath, 'r') as file:
    data_file = json.load(file)

# Plot graph
fig, ax = plt.subplots(2)
#fig.set_size_inches(12, 8)

# Make patches for legends
ch_patch = mpatches.Patch(color='red', label='Cheater alleles')
re_patch = mpatches.Patch(color='blue', label='Resistor alleles')

mt1_patch = mpatches.Patch(color='red', label='Type 1')
mt2_patch = mpatches.Patch(color='blue', label='Type 2')
#mt3_patch = mpatches.Patch(color='green', label='Type 3')

# Plot genotypes
ax[0].set_ylabel("Percentage of possible alleles")
ax[0].set_title("")
ax[0].set_ylim([0, 1])

ax[0].plot(data_file['x_axis_values'], data_file['mean_ch'], color = "red")
ax[0].fill_between(data_file['x_axis_values'], (np.array(data_file['mean_ch'])-np.array(data_file['ch_ci'])), (np.array(data_file['mean_ch'])+np.array(data_file['ch_ci'])), color='red', alpha=.1)
ax[0].plot(data_file['x_axis_values'], data_file['mean_res'], color = "blue")
ax[0].fill_between(data_file['x_axis_values'], (np.array(data_file['mean_res'])-np.array(data_file['res_ci'])), (np.array(data_file['mean_res'])+np.array(data_file['res_ci'])), color='blue', alpha=.1)
ax[0].legend(handles=[ch_patch, re_patch])

# Plot mating types
ax[1].set_ylabel("Ratio")
ax[1].set_title("Mating type distribution over time")
ax[1].set_ylim([0, 1])

ax[1].plot(data_file['x_axis_values'], data_file['mean_mt1'], color = "red")
ax[1].plot(data_file['x_axis_values'], data_file['mean_mt2'], color = "blue")
#ax[1].plot(data_file['x_axis_values'], data_file['mean_mt3'], color = "green")
ax[1].fill_between(data_file['x_axis_values'], (np.array(data_file['mean_mt1'])-np.array(data_file['mt1_ci'])), (np.array(data_file['mean_mt1'])+np.array(data_file['mt1_ci'])), color='red', alpha=.1)
ax[1].fill_between(data_file['x_axis_values'], (np.array(data_file['mean_mt2'])-np.array(data_file['mt2_ci'])), (np.array(data_file['mean_mt2'])+np.array(data_file['mt2_ci'])), color='blue', alpha=.1)
ax[1].fill_between(data_file['x_axis_values'], (np.array(data_file['mean_mt3'])-np.array(data_file['mt3_ci'])), (np.array(data_file['mean_mt3'])+np.array(data_file['mt3_ci'])), color='green', alpha=.1)
ax[1].legend(handles=[mt1_patch, mt2_patch])


ax[0].set_xlabel("Development cycles")
ax[1].set_xlabel("Development cycles")
plt.subplots_adjust(hspace=0.5)


# Plot lines indicating sexual cycles
if data_file["parameters"]["i_macrocyst"][0] != 0:
    for i in data_file["parameters"]["i_macrocyst"]:
        ax[0].axvline(x = i, color = "b", linestyle = "dashed", alpha = 0.05)
        ax[1].axvline(x = i, color = "b", linestyle = "dashed", alpha = 0.05)
        
else:
    for i in data_file["sex_cycle_list"]:
        ax[0].axvline(x = i, color = "b", linestyle = "dashed", alpha = 0.05)
        ax[1].axvline(x = i, color = "b", linestyle = "dashed", alpha = 0.05)
        
# Create canvas for gene pair graphs
gene_pairs = data_file["parameters"]["gene_pairs"]
fig, ax = plt.subplots(gene_pairs)

for x in range(0, gene_pairs):
    temp_str = "gene_pair" + str(x)

    # Set title and axis labels
    ax[x].set_title("Genotypes at gene pair " + str(x+1))
    ax[x].set_ylabel("# of individuals")
    ax[x].set_xlabel("Developmental cycles")

    # Plot mean gene pair specific data from .json data file - this shows the mean absolute number of individuals with each geneotype at a specific gene pair
    ax[x].set_ylim([0, 10000])
    ax[x].plot(data_file['x_axis_values'], data_file['average_tracker_dict'][temp_str][0], color = "red")
    ax[x].plot(data_file['x_axis_values'], data_file['average_tracker_dict'][temp_str][1], color = "blue")
    ax[x].plot(data_file['x_axis_values'], data_file['average_tracker_dict'][temp_str][2], color = "green")

    # Draw confidence interval data around mean data points
    ax[x].fill_between(data_file['x_axis_values'], (np.array(data_file['average_tracker_dict'][temp_str][0])-np.array(data_file['ci_tracker_dict'][temp_str][0])), (np.array(data_file['average_tracker_dict'][temp_str][0])+np.array(data_file['ci_tracker_dict'][temp_str][0])), color='red', alpha=.1)
    ax[x].fill_between(data_file['x_axis_values'], (np.array(data_file['average_tracker_dict'][temp_str][1])-np.array(data_file['ci_tracker_dict'][temp_str][1])), (np.array(data_file['average_tracker_dict'][temp_str][1])+np.array(data_file['ci_tracker_dict'][temp_str][1])), color='blue', alpha=.1)
    ax[x].fill_between(data_file['x_axis_values'], (np.array(data_file['average_tracker_dict'][temp_str][2])-np.array(data_file['ci_tracker_dict'][temp_str][2])), (np.array(data_file['average_tracker_dict'][temp_str][2])+np.array(data_file['ci_tracker_dict'][temp_str][2])), color='green', alpha=.1)

    # Draw sex cycle indicators
    if data_file["parameters"]["i_macrocyst"][0] != 0:
        for i in data_file["parameters"]["i_macrocyst"]:
            ax[x].axvline(x = i, color = "b", linestyle = "dashed", alpha = 0.05)
    else:
        for i in data_file["sex_cycle_list"]:
            ax[x].axvline(x = i, color = "b", linestyle = "dashed", alpha = 0.05)
plt.subplots_adjust(hspace=1.5)

plt.show()