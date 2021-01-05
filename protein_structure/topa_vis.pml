# SCRIPT NAME:  topa_vis.pml
# Dan Deatherage  June 2020
# Comments for script file
# Rare variant visualization of topA mutation locations
# To run the script in pymol:   Select File --> Run Script--> then browse to the topa_vis.pml file


# basic set up
reinitialize
fetch 1mw8, async=0
set seq_view, 1
bg_color white
hide nonbonded
util.cbc

# create selections for different molecule types
sele so4, resn SO4
sele tmp, resn TMP
sele DNA, chain Y
hide cartoon, DNA
show sticks, DNA
util.cbaw so4 + tmp
sele topA, Chain X &! (resn SO4 | resn TMP)
color marine, DNA
hide cartoon, tmp
show spheres, tmp
color black, tmp

# view for scenes
set_view (\
     0.953917921,   -0.120496064,    0.274816900,\
    -0.243878260,    0.222278744,    0.943990171,\
    -0.174833119,   -0.967508316,    0.182648405,\
     0.000000000,    0.000000000, -337.424804688,\
    44.090461731,   26.789920807,   30.533090591,\
   266.028289795,  408.821319580,  -20.000000000 )


scene topA_vs_DNA, store

# create selections for different domains
sele D1, resi 1-35 + resi 83-157
sele D2, resi 214-281 + resi 406-471
sele D3, resi 282-405
sele D4, resi 36-82 + resi 158-213 + resi 472-589
# following domains are beyond 1mw8 length
# sele D5 resi 590-NA
# sele D6 resi NA-700
# sele D7 resi 701-750
# sele D8 resi 751-NA
# sele D9 resi NA-865

set cartoon_transparency, 0.4, 1mw8

# color things in light colors (reserving Green for mutations)
util.cbaw topA
util.cbab D1
util.cbam D2
util.cbas D3
util.cbay D4

scene Domains, store


# select mutations (taken from mutations_for_pymol.py output)

## Original
#sele topA_Muts, resi 34 + resi 36 + resi 74 + resi 76 + resi 109 + resi 115 + resi 127 + resi 134 + resi 139 + resi 151 + resi 179 + resi 187 + resi 191 + resi 192 + resi 205 + resi 286 + resi 293 + resi 323 + resi 367 + resi 391 + resi 496 + resi 560 + resi 567 + resi 586
# below are additional mutations in regions not modeled by 1mw8
# + resi 617 + resi 622 + resi 761 + resi 774 + resi 778 + resi 780 + resi 793

## Updated for new mutations 2020-12-07
sele topA_Muts, resi 34 + resi 35 + resi 36 + resi 74 + resi 76 + resi 109 + resi 115 + resi 121 + resi 127 + resi 128 + resi 134 + resi 139 + resi 151 + resi 165 + resi 179 + resi 187 + resi 191 + resi 192 + resi 193 + resi 205 + resi 286 + resi 290 + resi 293 + resi 323 + resi 328 + resi 367 + resi 391 + resi 494 + resi 496 + resi 560 + resi 567 + resi 586
# below are additional mutations in regions not modeled by 1mw8
# + resi 617 + resi 622 + resi 761 + resi 774 + resi 778 + resi 780 + resi 793

util.cbag topA_Muts

set cartoon_transparency, 0.0, topA_Muts
scene mutations_cartoon, store


show spheres, topA_Muts
hide everything, so4
hide everything, tmp
scene mutations_spheres, store

# set of scenes to show each domain's mutations indvidually
hide (topA &! (D1 | tmp | so4 | DNA))
scene D1_mutations, store


scene mutations_spheres, recall
hide (topA &! (D2 | tmp | so4 | DNA))
scene D2_mutations, store


scene mutations_spheres, recall
hide (topA &! (D3 | tmp | so4 | DNA))
scene D3_mutations, store


scene mutations_spheres, recall
hide (topA &! (D4 | tmp | so4 | DNA))
scene D4_mutations, store


scene mutations_spheres, recall
hide (D1 | D2 | D3 | D4)
scene noDomain_mutations, store

scene mutations_spheres, recall
sele nearTMP, tmp around 5
hide (topA &! (nearTMP | tmp | so4 | DNA))
scene muts_wi_5A_tmp, store


scene mutations_spheres, recall
sele nearDNA, DNA around 5
hide (topA &! (nearDNA | tmp | so4 | DNA))
scene muts_wi_5A_DNA, store


# labeling commands from pykf file for reference. Determine what lables are appropriate.
# would be easier if knew how to move labels via command line

# # labeling. Note that the labels are moved by hand for .pse file
# # unknown how to modify positions of the labels in this script
# pseudoatom L1_domainB
# label L1_domainB, "Domain B"
# pseudoatom L2_domainC
# label L2_domainC, "Domain C"
# pseudoatom L3_domainA
# label L3_domainA, "Domain A"
# pseudoatom L4_monomer
# label L4_monomer, "Monomer"
# pseudoatom L5_C_ineractions
# label L5_C_ineractions, "Domain C 5A connections"
# pseudoatom L6_A_interactions
# label L6_A_interactions, "Domain A 5A connections"
# set seq_view, 0
# scene labeled, store
