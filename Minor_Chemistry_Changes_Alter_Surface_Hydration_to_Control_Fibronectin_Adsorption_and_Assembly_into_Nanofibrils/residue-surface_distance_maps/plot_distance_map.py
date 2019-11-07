#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Visualise the distance map representing residues-surface distances calculated using the
accompanying python script.

Please cite the following publication if you find this script useful:
Bieniek, M. K. et al. (2019) ‘Minor Chemistry Changes Alter Surface Hydration to Control Fibronectin Adsorption and Assembly into Nanofibrils’, Advanced Theory and Simulations. John Wiley & Sons, Ltd, p. 1900169. doi: 10.1002/adts.201900169.

Bibtex:

@article{Bieniek2019,
author = {Bieniek, Mateusz K. and Llopis‐Hernandez, Virginia and Douglas, Katie and Salmeron‐Sanchez, Manuel and Lorenz, Christian D.},
doi = {10.1002/adts.201900169},
file = {:home/dresio/Downloads/de/Bieniek{\_}et{\_}al-2019-Advanced{\_}Theory{\_}and{\_}Simulations.pdf:pdf},
issn = {2513-0390},
journal = {Advanced Theory and Simulations},
keywords = {fibrillogenesis,fibronectin,material‐driven fibrillogenesis,molecular dynamics,surface hydration},
month = {oct},
pages = {1900169},
publisher = {John Wiley {\&} Sons, Ltd},
title = {{Minor Chemistry Changes Alter Surface Hydration to Control Fibronectin Adsorption and Assembly into Nanofibrils}},
url = {https://onlinelibrary.wiley.com/doi/abs/10.1002/adts.201900169},
year = {2019}
}

Please create an issue on github if you need any assistance or have any suggestions.
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm


# paths which have to be adjusted by the user
data_filepath = "data_distance_map.dat"
plot_output_filepath = "distance_map.png"

# general configuration for matplotlib
matplotlib.rcParams.update({'font.size': 8})
plt.figure(figsize=(6, 2))
# just one plot
plt.subplot(1, 1, 1)
plt.title('EA10 Title')
plt.xlabel('Time (ns)')
plt.ylabel('Residues')

# read the residue names first
# this is format dependant and depends on how it was saved in the first place
resnames = open(data_filepath).readline().split('time(ps)')[1].split()
resnames = [n.split('.')[0] + ' ' + n.split('.')[1] for n in resnames]

# load the data with numpy
# plot only every 10th frame (every 10th row), in my case that's a step of 1 ns
data = np.loadtxt(data_filepath, comments='#')
assert len(resnames) == data.shape[1] - 1, \
    'The number of residues does not correspond to the number of data columns ' \
    '(no of residues + 1 column for time)'
# remove the time column
heatmap = data[:, list(range(1, len(resnames) + 1))]
# the x axis should be time, and the y axis should be residues
heatmap = heatmap.T
# plot while ignoring distances above 20 angstroms
# remember to adjust colour map: jet_r is jet is reversed
# list of colour maps: https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html
contactmap = plt.pcolormesh(heatmap, cmap=cm.jet_r, vmin=0, vmax=20)

# prepare places for the yticks (residue names)
tick_indices = np.array(list(range(0, len(resnames) + 1, 20)))
# shift each label to apply in the middle of its "square"
# adjust the fontsize of the y ticks
plt.yticks(tick_indices + 0.5, [resnames[i] for i in tick_indices], fontsize=6)

# create the legend bar in the figure, the ticks correspond to the distances
cbar = plt.colorbar(contactmap, ticks=[0, 5, 10, 15, 20], fraction=0.02, pad=0.04)
cbar.ax.set_yticklabels(['0$\\rm \AA$', '5$\\rm \AA$', '10$\\rm \AA$', '15$\\rm \AA$', '>20$\\rm \AA$'])

plt.tight_layout()
plt.savefig(plot_output_filepath, dpi=300)
# plt.show()

