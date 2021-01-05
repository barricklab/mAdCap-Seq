# adaptome-capture

Scripts and processed data files to accompany:

Deatherage, D.E., Barrick, J.E. **Profiling the initial burst of beneficial genetic diversity in clonal cell populations to anticipate evolution**. https://doi.org/10.1101/2020.07.10.196832

## Overview

`reference` contains sequence files for the targeted genes of interest
and the _E. coli_ REL606 genome with these regions masked out.
`consensus_read_generation` has the script for using molecular indexes
to perform error correction on raw Illumina reads. After running
[_breseq_](https://github.com/barricklab/breseq) to predict genetic variants using these files as input, the
output is converted into a format with the reads supporting the
reference versus variant alleles using `breseq_postprocessing`. Then
`trajectory_analysis` contains the main scripts for filtering and
analyzing the trajectories of mutation frequencies. `LTEE-compare` and
`protein_structure` contain scripts and information for further analyzing
the sets of predicted mutations.
