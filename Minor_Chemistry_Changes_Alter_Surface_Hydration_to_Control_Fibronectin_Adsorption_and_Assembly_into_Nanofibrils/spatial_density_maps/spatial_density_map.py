#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
In order to understand the environment composition / density with respect
to the reference molecule, superimpose a selected structure and reorient the
environment accordingly. Save the superimposed structures with the
surrounding water molecules in a .pdb file.

------- Post-Processing Instructions -------

After the superimposition and extraction of the surrounding molecules,
extract the water oxygens from the output file, using for example this command:
grep OW sdm.pdb > sdm_ow.pdb

Then, load the sdm.pdb file which includes all the water oxygens. Next use
 VMD->Analysis->Volmap to generate a sdm.dx map (selection=all).

Then, load a single PEAC chain, and the load the sdm.dx file. Change the
visualisation in Graphics->Representations for the sdm.dx to Isosurface
and use the Draw: Solid Surface. Modify the isovalue to visualise different
densities.

------- Post-Processing Instructions -------

Instructions: read carefully the script and update the configuration accordingly.

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

Thanks to:
MDAnalysis: https://www.mdanalysis.org/
This hydration analysis was inspired by ANGULA created by Luis Carlos Pardo:
https://gcm.upc.edu/en/members/luis-carlos/angula/ANGULA
"""

import MDAnalysis as mda
from MDAnalysis.analysis.align import rotation_matrix
import matplotlib.pyplot as plt


def fix_pbc(coord, pbc):
    """
    Move the coordinate to the right PBC.
    Assumes that the translation to the origin has been done.
    And that the selection is smaller than half the PBC.
    :param coord: coordinate along the pbc axis
    :param pbc: pbc size
    """
    if coord > pbc / 2:
        return coord - pbc
    elif coord < -(pbc / 2):
        return coord + pbc

    return coord

# load the simulation in MDAnalysis
u = mda.Universe('sam_hydrated.gro', 'sam_hydrated_step10ns.xtc')
# the example input file is a Self-assembled Monolayer (SAM) with 252 residues namd PEAC

# use the first frame of our trajectory as the reference structure
# which will be used to superimpose other frames
# select four atoms which include an ester as our reference
ea_template = u.select_atoms('resname PEAC and name C19 O2 O1 C20 and resid 1', updating=False)
# shift the origin to be the position of the C20 atom
ea_template = ea_template.atoms.translate(-ea_template.select_atoms('name C20')[0].position)
ea_ref_pos = ea_template.positions

# define the output .pdb structure where the superimposed frames will be stored. Note: it has to be a .pdb
# due to the flexible number of water molecules found.
output = 'sdm.pdb'

with mda.Writer(output) as W:
    rmsds = []
    # for each frame
    for ts in u.trajectory:
        print("Time (ns):", ts.time / 1000)
        # for each molecule in the self-assembled monolayer
        for res in u.select_atoms('resname PEAC').residues:
            # select the residue together with the water oxygens within 5.8 of the C20 atom in the residue
            sel = u.select_atoms('resid %d or (name OW and around 5.8 (resid %d and name C20))'
                                 % (res.resid, res.resid)).residues.atoms
            original_positions = sel.positions
            structural_template = sel.select_atoms('name C19 O2 O1 C20')

            # set the position of c20 to be the origin
            c20 = structural_template.select_atoms('name C20')[0]
            sel.translate(-c20.position)

            # correct the PBC
            for atom in sel:
                corrected_x = fix_pbc(atom.position[0], u.dimensions[0])
                corrected_y = fix_pbc(atom.position[1], u.dimensions[1])
                corrected_z = fix_pbc(atom.position[2], u.dimensions[2])

                atom.position = (corrected_x, corrected_y, corrected_z)

            # find the rotation matrix necessary to superimpose the structure against the reference
            R, rmsd = rotation_matrix(structural_template.atoms.positions, ea_ref_pos)
            # monitor the quality of the superimposition
            rmsds.append(rmsd)
            # apply the rotation matrix to the molecule, including the surrounding water
            sel.rotate(R)
            # save the superimposed atoms in a .pdb format
            W.write(sel)


plt.title('Superimposed against reference')
plt.ylabel('RMSD')
plt.hist(rmsds)
plt.show()