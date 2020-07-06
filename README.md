# adaptome-capture

Scripts and processed data files to accompany:

Deatherage, D.E., Barrick, J.E. **Profiling the initial burst of
beneficial genetic diversity to anticipate evolution of a cell
population**

## Overview

`reference` contains sequence files for the targeted genes of interest
and the _E. coli_ REL606 genome with these regions masked out.
`consensus_read_generation` has the script for using molecular indexes
to perform error correction on raw Illumina reads. After running
_breseq_ to predict genetic variants using these files as input, the
output is converted into a format with the reads supporting the
reference versus variant alleles using `breseq_output_conversion`. Then
`trajectory_analysis` contains the main scripts for filtering and
analyzing the trajectories of mutation frequencies. `LTEE-compare` and
`protein_structure` contain scripts for the final analyses of the sets
of predicted mutations.
