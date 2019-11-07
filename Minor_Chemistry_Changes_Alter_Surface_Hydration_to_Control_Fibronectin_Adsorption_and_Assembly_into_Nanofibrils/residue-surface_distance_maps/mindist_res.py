#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculate the minimum distance from the heavy atoms of each residue and the surface.
Save the results in a file with a tabular format with the following format:
# time(ps) res1aaa.res1id.res1index res2aaa.res2id.res2index

Instructions: read carefuly the script and update the configuration accordingly.

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

import MDAnalysis
from MDAnalysis.analysis.distances import distance_array
import numpy as np

# open the trajectory with MDAnalysis
u = MDAnalysis.Universe('500ns_protCent_pbcMol.gro',
                        '500ns_centCA_pbcMol_step100ps.xtc')
# select the atoms in the substrate that are found on the surface
substrate = u.select_atoms("resname EA and (name C11 O1 O2 C12 C13)")
print ('Selected', substrate)
# the name of the file where the data will be saved
output_filename = 'data_distance_map.dat'

# select the protein
protein = u.select_atoms('protein')

# create column titles with the format "time, res1, res2, res3"
column_titles = ['time(ps)', ]
# preselect the heavy atoms for each residue (optimisation)
residues = []
for i, residue in enumerate(protein.residues, start=1):
    column_titles.append(residue.resname + '.' + str(residue.resid) + '.' + str(i))
    residues.append(residue.atoms.select_atoms('not name H*'))

assert len(protein.residues) == len(residues), 'There should be heavy atom indices for each residue'

data = []
# adjust the step of the trajectory
for ts in u.trajectory[::1]:
    row = [ts.time, ]
    for res_hatoms in residues:
        # find the minimum distance from the heavy residues to the substrate
        mindst = np.min(distance_array(res_hatoms.positions, substrate.positions, box=ts.dimensions))
        row.append(mindst)
    data.append(row)
    print ("Done timeframe (ns):", ts.time / 1000)

# save the data with numpy
np.savetxt(output_filename, data, fmt='%.2f', header=' '.join(column_titles))
