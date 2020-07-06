#!/usr/bin/env python3

"""
Small python file to:
 read in significant mutations used in figure generation
 identify mutations that would be useful to visualize in pymol (SNPs)
 convert from DNA position to AA position (note varying length flanking sequence)
 generate pymol selection command
"""

import os
import re

def ceiling_division(n, d):
    """
    math.ceil can have problems with 'large' integers as it converts from ints to
    floats and back. This isn't very fast and it may have rounding issues.

    This function adapts floor division and semi non-intuitive negative integer
    interaction wherein -1.1 floor division down to -2 not -1 as it is not front-
    end rounding but rather rounding down to next lowest integer.

    Used in converting DNA -> AA position function
    """

    assert n >= 0 and d >= 0, "ceiling_division function requires both numbers to be positive. in this case n = %s and d = %s causing assertion error" % (n,d)

    return -(-n // d)

def aa_position(dna_position=0, testing=False):
    """function to turn DNA based positioning into AA positioning"""
    if testing:
        print("DNA", "AA", sep="\t")
        for _ in range(10):
            print(_, ceiling_division(_, 3), sep="\t")
        print("note that 0 // 0 gives 0 and hence would not make sense in AA space")

    assert dna_position != 0, "DNA position given as 0bp this will not convert to an AA. rerun aa_position in testing mode for more info."

    return ceiling_division(dna_position, 3)

dna_positions = []
for file_in in os.listdir("../Rare_variant_figure_data"):
    if file_in.endswith("significant_mutations.csv"):
        with open("../Rare_variant_figure_data/" + file_in, "r") as fin:
            for line in fin:
                if re.search("pykF.*nonsynonymous", line):
                    line=line.rstrip().split(",")
                    try:
                        pos = int(line[3]) - 1630  # as per 1400 bp flanking file but want bp 1 to be first base of cds
                    except(ValueError):  # in window file, position in column E not D
                        pos = int(line[4]) - 1630  # as per 1400 bp flanking file but want bp 1 to be first base of cds
                    assert 1 <= pos <= 1412, "Converted DNA position not within bounds of gene length %s" % pos
                    dna_positions.append(pos)
print(len(dna_positions), "total pykF mutations")
print(sorted(set(dna_positions)))
print(len(sorted(set(dna_positions))), "pykF mutations at unique DNA positions")

aa_positions = [aa_position(x) for x in dna_positions]
print(sorted(set(aa_positions)))
print(len(sorted(set(aa_positions))), "different pykF AA mutated")

print("sele pykF_Muts, resi " + " + resi ".join(str(x) for x in sorted(set(aa_positions))))


dna_positions = []
for file_in in os.listdir("../Rare_variant_figure_data"):
    if file_in.endswith("significant_mutations.csv"):
        with open("../Rare_variant_figure_data/" + file_in, "r") as fin:
            for line in fin:
                if re.search("topA.*nonsynonymous", line):
                    line=line.rstrip().split(",")
                    try:
                        pos = int(line[3]) - 1658  # as per 1400 bp flanking file but want bp 1 to be first base of cds
                    except(ValueError):  # in window file, postition in column E not D
                        pos = int(line[4]) - 1658  # as per 1400 bp flanking file but want bp 1 to be first base of cds
                    assert 1 <= pos <= 2597, "Converted DNA position not within bounds of gene length %s" % pos
                    dna_positions.append(pos)
print(len(dna_positions), "total topA mutations")
print(sorted(set(dna_positions)))
print(len(sorted(set(dna_positions))), "topA mutations at unique DNA positions")

aa_positions = [aa_position(x) for x in dna_positions]
print(sorted(set(aa_positions)))
print(len(sorted(set(aa_positions))), "different topA AA mutated")

print("sele topA_Muts, resi " + " + resi ".join(str(x) for x in sorted(set(aa_positions))))