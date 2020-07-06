# SCRIPT NAME:  pykf_vis.pml
# Dan Deatherage  May 2020
# Comments for script file
# Rare variant visualization of pykF mutation locations
# To run the script in pymol:   Select File --> Run Script--> then browse to the scritpname.pml file


# two structures seem different
# fetch 4yng, async=0
# fetch 1pky, async=0
# super 4yng, 1pky
# matching score 7737.786; Executive RMSD 41.403
# 4yng worked with going forward


# basic set up and emphasis for SO4 molecules
reinitialize
fetch 4yng, async=0
set seq_view, 1
bg_color white
hide nonbonded
util.cbc
sele so4, resn SO4
util.cbay so4

# View appropriate only if pseudo atoms are not added; they change size of seq_view window
# set_view (\
    -0.104747556,    0.792604387,    0.600672305,\
     0.164769471,   -0.581816316,    0.796455145,\
     0.980754495,    0.182399184,   -0.069653004,\
     0.000000000,    0.000000000, -389.156097412,\
   -32.694007874,  -13.068845749,   62.495571136,\
   230.159210205,  548.153015137,  -20.000000000 )

# Alt view appropriate if psuedo atoms were/are-to-be added  
set_view (\
    -0.103447556,    0.789930403,    0.604408443,\
     0.163329959,   -0.585928917,    0.793732107,\
     0.981133163,    0.180827633,   -0.068406247,\
     0.000000000,    0.000000000, -608.419799805,\
   -32.694007874,  -13.068845749,   62.495571136,\
   449.422943115,  767.416687012,  -20.000000000 )

scene two_tetramers_color_by_chain, store


set cartoon_transparency, 0.4, 4yng

# set up selections for different domains 
sele domainA, resi 1-70 + resi 171-345
sele domainB, resi 71-170
sele domainC, resi 352-470

# color each domain in light colors
util.cbaw 4yng
util.cbas domainA
util.cbay domainB
util.cbab domainC
scene color_by_domain, store


# select mutations (taken from mutations_for_pymol.py output)
sele pykF_Muts, resi 9 + resi 69 + resi 70 + resi 92 + resi 101 + resi 144 + resi 177 + resi 190 + resi 193 + resi 224 + resi 239 + resi 240 + resi 243 + resi 245 + resi 248 + resi 254 + resi 257 + resi 261 + resi 264 + resi 278 + resi 280 + resi 296 + resi 301 + resi 306 + resi 312 + resi 313 + resi 327 + resi 328 + resi 353 + resi 362 + resi 380 + resi 381 + resi 385 + resi 400 + resi 460
util.cbag pykF_Mut
set cartoon_transparency, 0.0, pykF_Muts
scene mutations_cartoon, store


show spheres, pykF_Muts
scene mutations_spheres, store

# begin to isolate different parts of the protein for viewing/labeling
hide (chain A &! (domainB | resn SO4))
hide (chain B &! (domainC | resn SO4))
hide (chain D &! (domainA | resn SO4))
hide (chain C &! resn SO4)
scene tetramer_vs_domain_mutations, store

# Currently set as atoms in adjacent monomers, NOT both sides of interaction
sele chainE_Contacts, chain E around 5
hide (chain F &! (chainE_Contacts | resn SO4))
hide (chain G &! (chainE_Contacts | resn SO4))
scene monomer_vs_momomerConnections_vs_domain_mutations, store

# i think could be adjusted by adding command like following 2 lines to make interface
# sele chain_interface, chainE_Contacts around 5 & chain E
#﻿hide (chain E &! chain_interface)


hide (so4)
util.cbc
scene mutations_by_chain_and_domain, store

# see comment on line 85? about doing as interface rather than interaction
hide (chain E)
scene mutations_by_domain_or_class, store

# recolor for next scene
util.cbas domainA
util.cbay domainB
util.cbab domainC

scene mutations_by_domain, store


# labeling. Note that the labels are moved by hand for .pse file
# unknown how to modify positions of the labels in this script
pseudoatom L1_domainB
label L1_domainB, "Domain B"
pseudoatom L2_domainC
label L2_domainC, "Domain C"
pseudoatom L3_domainA
label L3_domainA, "Domain A"
pseudoatom L4_monomer
label L4_monomer, "Monomer"
pseudoatom L5_C_ineractions
label L5_C_ineractions, "Domain C 5A connections"
pseudoatom L6_A_interactions
label L6_A_interactions, "Domain A 5A connections"
set seq_view, 0
scene labeled, store
