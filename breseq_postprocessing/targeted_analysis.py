#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 10:14:39 2013
options:
-m = all mutations
-m -e = all mutations and all evidence
-se = supporting evidence only
-m -ese = equivalent to -m
-m -ese -e = mutations and orphan evidence
@author: ded
V2 = added option for generational time points being added via additional file
V4 = cleanup for paper
"""

import os
import re
import argparse
import numpy
import sys


class FullPaths(argparse.Action):
    """Expand user- and relative-paths
    taken from https://gist.github.com/brantfaircloth/1252339 to allow relative path in argparse for gd directory identification
    """
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


parser = argparse.ArgumentParser(description="read in fastq file names, generate .gd files")
parser.add_argument("-g", "--gd", help="directory containing gd files to compare", action=FullPaths)
parser.add_argument("-l", "--log", help="log output filename")
parser.add_argument("-p", "--prefix", help="prefix for output names")
parser.add_argument("-m", "--mutations", action='store_true', default=False, help="evaluate all 3 letter mutation differences")
parser.add_argument("-e", "--evidence", action='store_true', default=False, help="evaluate all 2 letter mutation differences (no handling for MC and UN)")
parser.add_argument("-mob", "--mobile_element", action='store_true', default=False, help="evaluate mobile element mutations")
parser.add_argument("-del", "--deletion", action='store_true', default=False, help="evaluate deletions mutations")
parser.add_argument("-snp", "--single_nucleotide_polymorphism", action='store_true', default=False, help="evaluate single nucleotide polymorphisms mutations")
parser.add_argument("-ins", "--insertion", action='store_true', default=False, help="evaluate insertion mutations")
parser.add_argument("-amp", "--amplification", action='store_true', default=False, help="evaluate amplification mutations")
parser.add_argument("-sub", "--substitution", action='store_true', default=False, help="evaluate substitution mutations (untested)")
parser.add_argument("-inv", "--inversion", action='store_true', default=False, help="evaluate inversion mutations (untested)")
parser.add_argument("-con", "--conversion", action='store_true', default=False, help="evaluate conversion mutations(untested)")
parser.add_argument("-jc", action='store_true', default=False, help="evaluate JC evidence differences")
parser.add_argument("-ra", action='store_true', default=False, help="evaluate RA evidence differences")
parser.add_argument("-sd", "--strand_discrepancy", type=float, default=False, help="attempt to determine strand discrepancies of RA evidence with strand bias in excess of given value")
parser.add_argument("-se", "--supporting_evidence", action='store_true', default=False, help="include evidence used to support mutations in output (using this option will cause 2 JC evidence lines to be collected for each MOB) Maybe required when looking for any evidence of a mutation across other lines.")
parser.add_argument("-ese", "--exclude_supporting_evidence", action="store_true", default=False, help="do not include evidence that was used to call mutations. Assumed use in conjuncture with -m and -e")
parser.add_argument('gdnames', nargs='*', help="file names of .gd files to compare. CAN NOT BE USED WITH -g option")
parser.add_argument("-v", "--verbose", action='store_true', help="print to screen")
parser.add_argument("-ac", "--autocorrelate", action='store_true', help="gd files are of a time course or other autocorrelative process, use this information to determine additional false positives")
parser.add_argument("-acm", "--auto_correlate_most", default=0.6, type=float, help="minimum autocorrelation value for mutation/evidence to be considered valid if AT LEAST 50%% of samples have mutation/evidence")
parser.add_argument("-acs", "--auto_correlate_some", default=0.2, type=float, help="minimum autocorrelation value for mutation/evidence to be considered valid if LESS THAN 50%% of samples have mutation/evidence")
parser.add_argument("-t", "--time", default=None, help="tsv with sample names and generational times")

# show help if no arguments passed
if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()


_nsre = re.compile('([0-9]+)')  # taken from stack overflow to "natural sort"


def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(_nsre, s)]


def printing(stuff_to_print):
    """Function for printing to screen if verbose flag supplied, or log file if log file supplied"""
    if args.verbose:
        print stuff_to_print
    if args.log is not None:
        print>>open(args.log, "a"), stuff_to_print


# assertions to make sure correct options or option sets are selected, and any relevant warnings
assert not (args.gd is None and args.gdnames == []), "gd files or directory containing gd files must be specified"
assert args.gd is None or args.gdnames == [], "gd files or directory must be specified, not both"
assert args.prefix is not None, "prefix must be given to name output files using the -p flag"
assert not (args.verbose is False and args.output is None), "-v or -o FILENAME required to print logs"
if args.evidence:
    printing("!!!WARNING!!!\nMC and UN evidence have no frequency associated with them, and is not displayed.")

# generate files_to_compare list of gd files to be compared regardless of how they were entered.
if args.gd:
    os.chdir(args.gd)  # TODO learn how to list directories (relative and absolute) without changing to that directory
    files_to_compare = [x for x in os.listdir(".") if x.endswith(".gd")]
else:
    files_to_compare = [x for x in args.gdnames if x.endswith(".gd")]
    assert len(files_to_compare) == len(args.gdnames), "Non Genome Diff file provided: %s" % ([x for x in args.gdnames if not x.endswith(".gd")])

# generate list of mutations to examine:
mutation_analysis_list = []
evidence_analysis_list = []
if args.mobile_element:
    mutation_analysis_list.append("MOB")
if args.deletion:
    mutation_analysis_list.append("DEL")
if args.single_nucleotide_polymorphism:
    mutation_analysis_list.append("SNP")
if args.insertion:
    mutation_analysis_list.append("INS")
if args.amplification:
    mutation_analysis_list.append("AMP")
if args.substitution:
    mutation_analysis_list.append("SUB")
if args.inversion:
    mutation_analysis_list.append("INV")
if args.conversion:
    mutation_analysis_list.append("CON")
if args.jc:
    evidence_analysis_list.append("JC")
if args.ra:
    evidence_analysis_list.append("RA")
assert not (mutation_analysis_list is not [] and args.supporting_evidence), "selected individual mutations to analyze in supporting evidence mode, this is not coded, justify the need and can be done"

# variable initiation for appends and trys
master_mutation_dict = {}
master_mutation_list = []

# generate list of mutations to exclude from analysis
# TODO soft code this based on reading in .gd file and using only locations, not names
excluded_mutations_list = []
excluded_mutations_list.append("REL606.JC.1.1.REL606.4629812.n1.REF.2")  # circularized genome
excluded_mutations_list.append("REL606.JC.1327760.n1.REL606_topa.1.1.REF.2")  # topA start
excluded_mutations_list.append("REL606.JC.1333417.1.REL606_topa.5655.n1.REF.2")  # topA end
excluded_mutations_list.append("REL606.JC.1731333.n1.REL606_pykf.1.1.REF.2")  # pykF start
excluded_mutations_list.append("REL606.JC.1735777.1.REL606_pykf.4442.n1.REF.2")  # pykF end
excluded_mutations_list.append("REL606.JC.3759033.n1.REL606_spot.1.1.REF.2")  # spoT start
excluded_mutations_list.append("REL606.JC.3764265.1.REL606_spot.5230.n1.REF.2")  # spoT end
excluded_mutations_list.append("REL606.JC.4098498.n1.REL606_hslu.1.1.REF.2")  # hslU start
excluded_mutations_list.append("REL606.JC.4103237.1.REL606_hslu.4737.n1.REF.2")  # hslU end
excluded_mutations_list.append("REL606.JC.4139250.n1.REL606_yijc.1.1.REF.2")  # yijC start
excluded_mutations_list.append("REL606.JC.4143032.1.REL606_yijc.3780.n1.REF.2")  # yijC end
excluded_mutations_list.append("REL606.JC.4200334.n1.REL606_iclr.1.1.REF.2")  # iclR start
excluded_mutations_list.append("REL606.JC.4204045.1.REL606_iclr.3709.n1.REF.2")  # iclR end
excluded_mutations_list.append("REL606.JC.4614128.n1.REL606_nadr.1.1.REF.2")  # nadR start
excluded_mutations_list.append("REL606.JC.4618161.1.REL606_nadr.4031.n1.REF.2")  # nadR end
excluded_mutations_list.append("REL606.JC.472228.n1.REL606_ybal.1.1.REF.2")  # yabL start
excluded_mutations_list.append("REL606.JC.476762.1.REL606_ybal.4532.n1.REF.2")  # ybaL end

if args.time is not None:
    generation_dictionary = {}
    with open(args.time, "r") as F:
        for line in F:
            line = line.rstrip().split('\t')
            assert line[0] not in generation_dictionary, "duplicate sample listing in generation file %s" % line
            assert line[1] not in generation_dictionary.values(), "duplicate generational listing in generation file %s" % line
            generation_dictionary[line[0]] = float(line[1])

# generate master dictionaries
for files in files_to_compare:
    assert files.endswith(".gd"), "non-gd file found. This should never trigger as previous code should have eliminated this file: %s" % files
    mutation_supporting_evidence = []
    with open(files, "r") as F:
        if args.time is None:
            generation = files.replace(".gd", "generation")
        else:
            generation = generation_dictionary[files.replace(".gd", "")]
            master_mutation_dict[generation] = {}
        for line in F:
            line = line.rstrip().split('\t')
            if re.match("^#", line[0]) or line[0] == '':
                if re.match("^#=TIME", line[0]):
                    assert args.time is None, "sample and generational time file provided for .gd file containing generational time information"
                    assert re.search("generation", generation), "multiple #=TIME lines found.\n%s\n%s" % (generation, line)
                    assert generation not in master_mutation_dict, "Multiple files from same time point given"
                    generation = float(line[1])
                    master_mutation_dict[generation] = {}
                continue  # ignore other meta data, and blank lines check next line
            # elif re.search("rejected=", str(line)):  # turned off after conversation with Jeff, will deal with all post filtering later
                # continue  # evidence was rejected for various reason, therefore evidence should be ignored, regardless of other considerations
            elif line[0] in ["UN", "MC"]:
                continue  # UN and MC hav no frequencies or coverage associated with them, and currently no functionality for overlapping regions. therefore skip
            elif line[0] in mutation_analysis_list or (re.match("^[A-Z]{3}$", line[0]) and (args.mutations or args.supporting_evidence)):  # must come before 2 letter codes, and assumes standard .gd output of evidence after mutations
                mutation_supporting_evidence.extend([int(item) for item in line[2].split(",")])
                if args.supporting_evidence:  # if supporting evidence on, only supporting evidence displayed  TODO is this what you want?
                    continue
            elif int(line[1]) in mutation_supporting_evidence and args.exclude_supporting_evidence:  # evidence is supporting mutation, and exclude supporting evidence was selected
                continue
            elif int(line[1]) in mutation_supporting_evidence and args.supporting_evidence:  # this is required to avoid the "unanalyzed catch" at the end
                pass
            elif line[0] in evidence_analysis_list or (re.match("^[A-Z]{2}$", line[0]) and args.evidence):  # specific evidence listed to be analyzed, or all evidence to be analyzed
                pass
            else:  # mutation or evidence not selected for analysis
                continue
            line[3] = line[3].replace("-", "_")  # reference may contain hyphen which R converts to . and messes up delimination
            if line[0] == "SNP":
                mutation_name = line[3] + ".SNP." + line[4] + "." + line[5]  # reference.SNP.pos.new_base
            elif line[0] == "DEL":
                mutation_name = line[3] + ".DEL." + line[4] + "." + line[5]  # reference.DEL.pos.length
            elif line[0] == "AMP":
                mutation_name = line[3] + ".AMP." + line[4] + "." + line[5] + "." + line[6]  # reference.AMP.pos.length.copy_number
            elif line[0] == "INS":
                mutation_name = line[3] + ".INS." + line[4] + "." + line[5] + "." + str(line).split("insert_position=")[1].split("'")[0]  # reference.INS.pos.new_base.position
            elif line[0] == "MOB":
                mutation_name = line[3] + ".MOB." + line[4] + "." + line[5] + "." + line[6] + "." + line[7]  # reference.MOB.position.IS_type.orientation.repeated_base_count     ## may run into problems with ins_start= etc
            elif line[0] == "RA":
                if line[7] == ".":
                    line[7] = "_"  # "_" used rather than source "." as R is using . as delimiter This corresponds to a deletion RA
                if line[6] == ".":
                    line[6] = "_"  # "_" used rather than source "." as R is using . as delimiter This corresponds to an insertion RA
                mutation_name = line[3] + ".RA." + line[4] + "." + line[6] + "." + line[7] + "." + line[5]  # reference.RA.position.old_base.new_base.new_position
            elif line[0] == "JC":
                if line[5] == "-1":  # R converts - to . therefore change - to n to denote negative orientation
                    line[5] = "n1"
                if line[8] == "-1":  # R converts - to . therefore change - to n to denote negative orientation
                    line[8] = "n1"
                line[6] = line[6].replace("-", "_")  # reference may contain hyphen which R converts to . and messes up delimination
                mutation_name = line[3] + ".JC." + line[4] + "." + line[5] + "." + line[6] + "." + line[7] + "." + line[8] + ".REF."  # reference1.JC.position.orientation.reference2.position2.orientation2.REF.number_of_reference_counts
                refs_to_count = 2  # a maximum of 2 references can be counted
                if "side_1_read_count=NA" in line:
                    refs_to_count -= 1
                if "side_2_read_count=NA" in line:
                    refs_to_count -= 1
                mutation_name += str(refs_to_count)
            else:
                assert False, "unknown mutation/evidence type attempting to be handled\n%s" % line
            if re.search("gene_name=.*snp_type=", str(line)):
                for thing in line:
                    if re.match("^gene_name=", thing) or re.match("^snp_type=", thing):
                        mutation_name += "." + thing
            elif re.search("gene_name=", str(line)):
                for thing in line:
                    if re.match("^gene_name=", thing):
                        mutation_name += "." + thing + "." + "snp_type=indel"
                        break
            # elif re.search("deleted=1", str(line)):
            #     assert re.search("RA.*REL606_pykf', '1919", str(line)), "expected deleted=1 only to be found in specific RA example when gene_name missing.\n%s\t%s" % (line, files)
            #     mutation_name += ".gene_name=pykF.snp_type=indel"
            else:
                if "RA" in line[0]:
                    mutation_name += ".unknown_gene.snp_type=assumed_SNP"
                #assert "RA" not in line[0], "gene_name field not found in RA line.\n%s\n%s" % (line, files)  # current assumption is that RA lines have to have a gene_name field and an optional snp_type field, those lines lacking a snp_type field should be indels  #not true in plasmid seq
            # following if else block existed to address a breseq bug regarding junction overlap lengths, believed to be corrected as of 2/13/15 (email correspondance) therefore only else statement to be executed
            # if mutation_name == "REL606.JC.3894996.n1.REL606_nadr.1546.1.REF.1" and files == "DED70.gd" and mutation_name in master_mutation_dict[generation]:
            #    continue
            # else:
            assert mutation_name not in master_mutation_dict[generation], "duplicated mutation name has been generated for a single sample.\n%s\t%s" % (mutation_name, files)
            # if mutation_name in master_mutation_dict[generation]:
            #    print "%s\t%s" % (mutation_name, files)
            if mutation_name in excluded_mutations_list:
                continue
            try:
                master_mutation_dict[generation][mutation_name] = {"var_cov": "NA", "total_cov": "NA", "frequency": float(str(line).split("'frequency=")[1].split("'")[0])}  # default to NA for coverage (ie for mutations), store frequency
            except ValueError:
                assert ("frequency=NA" in line and "side_1_read_count=NA" in line and "side_2_read_count=NA") or ("frequency=NA" in line and "side_1_read_count=0" in line and "side_2_read_count=0"), "frequency could not be converted to a float, assumed that frequency must be listed as NA caused by both side1 and side2 read counts being 0 or NA.\n%s\n%s" % (line, files)
                master_mutation_dict[generation][mutation_name] = {"var_cov": "NA", "total_cov": "NA", "frequency": "NA"}  # default to NA for coverage (ie for mutations), Frequency listed as NA and can't be converted
            master_mutation_list.append(mutation_name)
            if line[0] == "RA":
                print line
                master_mutation_dict[generation][mutation_name]["var_cov"] = sum([int(x) for x in str(line).split("new_cov=")[1].split("'")[0].split("/")])
                master_mutation_dict[generation][mutation_name]["total_cov"] = sum([int(x) for x in str(line).split("total_cov=")[1].split("'")[0].split("/")])
            if line[0] == "JC":
                master_mutation_dict[generation][mutation_name]["var_cov"] = int(str(line).split("new_junction_read_count=")[1].split("'")[0])
                master_mutation_dict[generation][mutation_name]["total_cov"] = int(str(line).split("new_junction_read_count=")[1].split("'")[0])
                if "side_1_read_count=NA" not in line:
                    master_mutation_dict[generation][mutation_name]["total_cov"] += int(str(line).split("side_1_read_count=")[1].split("'")[0])
                if "side_2_read_count=NA" not in line:
                    master_mutation_dict[generation][mutation_name]["total_cov"] += int(str(line).split("side_2_read_count=")[1].split("'")[0])
                # following block based on assumption, and troubleshooting assumption that if 0 reads (NA) were detected for a side, it was because the region was redundant. 100% junctions also cause the same effect
                # if side1 == "NA":
                #     # assert "side_1_redundant=1" in line, "%s\n%s" % (line,files)
                #     if "side_1_redundant=1" not in line:
                #         print "%s\t%s" % (line, files)
                # else:
                #     total.append(int(side1))
                # if side2 == "NA":
                #     assert "side_2_redundant=1" in line
                # else:
                #     total.append(int(side2))
                # try:
                #     total = (sum(total) + master_mutation_dict[generation][mutation_name]["var_cov"]) / len(total)
                # except ZeroDivisionError:
                #     print files
                #     for thing in  line:
                #         print thing
                #     assert False
                # master_mutation_dict[generation][mutation_name]["total_cov"] = total

master_mutation_list = sorted(list(set(master_mutation_list)))  # generate non-redundant list of mutations so that each column in final output will be unique.. sorting allows grouping based on references, this may fall apart on junctions as potentially from different references
generation_list = sorted(master_mutation_dict.keys())  # intended to generate a list of generations used in analysis to facilitate autocorrelation analysis

if args.autocorrelate:
    printing("12-4-14, autocorrelation falling out of favor at this stage of pipeline in favor of using it or other filters later in process")
    final_mutation_dict = dict((gen, {}) for gen in generation_list)  # populate new dictionary with each existing generation and and empty dictionary
    for mutation in master_mutation_list:  # each mutation is separate
        valid_mut = False  # assume mutation will not autocorrelate
        # generate list of autocorrelation frequencies
        to_auto_correlate = []
        for i in xrange(1, len(generation_list)):  # can this be folded into next for loop?
            assert float(generation_list[i]) > float(generation_list[i - 1]), "Generation list not in increasing order, autocorrelation will not work.\n%s" % generation_list
        for generation in generation_list:
            try:
                to_auto_correlate.append(float(master_mutation_dict[generation][mutation]['frequency']))  # if mutation is present at given generation append that frequency to the correlation list
            except KeyError:
                to_auto_correlate.append(float(0))  # if mutation is not present at given generation, append 0 to the correlation list
            except ValueError:
                assert master_mutation_dict[generation][mutation]['frequency'] is "NA"
                to_auto_correlate.append(float(0))  # if mutation frequency is not present, it is treated as 0 ... this may not be correct
        t = 0.0
        b = 0.0
        for i in xrange(0, (len(to_auto_correlate) - 1)):
            t += (to_auto_correlate[i] - numpy.mean(to_auto_correlate)) * (to_auto_correlate[i + 1] - numpy.mean(to_auto_correlate))
            b += (to_auto_correlate[i] - numpy.mean(to_auto_correlate)) ** 2
        b += (to_auto_correlate[-1] - numpy.mean(to_auto_correlate)) ** 2
        # block is old method of autocorrelation, don't believe it is relevant any longer assertion added to verify truth
        assert sum(value < 0 for value in to_auto_correlate) is 0, "negative frequency detected.\n%s" % to_auto_correlate
        #     if re.search("NA", str(to_auto_correlate[i])):  # not sure that "NA" can currently be entered into the frequency key value
        #         to_auto_correlate[i] = 0
        #     t += (abs(float(to_auto_correlate[i])) - numpy.mean(to_auto_correlate)) * (abs(float(to_auto_correlate[i + 1])) - numpy.mean(to_auto_correlate))
        #     b += (abs(float(to_auto_correlate[i])) - numpy.mean(to_auto_correlate)) ** 2
        # b += (abs(float(to_auto_correlate[-1])) - numpy.mean(to_auto_correlate)) ** 2
        auto_correlation_coefficient = t / b
        if sum(value > 0 for value in to_auto_correlate) / float(len(to_auto_correlate)) >= 0.5:  # more than 50% of samples are mutated
            if auto_correlation_coefficient >= args.auto_correlate_most:
                valid_mut = True
        else:  # mutation occurs in less than 50% of samples
            if auto_correlation_coefficient >= args.auto_correlate_some:
                valid_mut = True
        if valid_mut:
            for generation in generation_list:
                assert mutation not in final_mutation_dict[generation], 'duplicate mutations generated?'  # should never catch as same assertion present in master_mutation_dictionary
                try:
                    final_mutation_dict[generation][mutation] = master_mutation_dict[generation][mutation]  # TODO check if this is circular reference, or if its a new instance
                except KeyError:  # mutation not found in given generation originally
                    pass
else:
    final_mutation_dict = master_mutation_dict  # this makes final_mutation_dict reference existing master_mutation_dict if autocorrelation not happening

with open(args.prefix + ".variant.tsv", "w") as Var, open(args.prefix + ".total.tsv", "w") as Tot, open(args.prefix + ".freq.tsv", "w") as Freq:
    print>>Var, "\t".join(map(str, ["generation"] + master_mutation_list))  # print header line
    print>>Tot, "\t".join(map(str, ["generation"] + master_mutation_list))  # print header line
    print>>Freq, "\t".join(map(str, ["generation"] + master_mutation_list))  # print header line
    for generation in generation_list:
        var_to_print = [generation]
        tot_to_print = [generation]
        freq_to_print = [generation]
        for mutation in master_mutation_list:
            try:
                var_to_print.append(final_mutation_dict[generation][mutation]['var_cov'])
            except KeyError:  # KeyError should only be raised when mutation not present at generation
                var_to_print.append("NA")
            try:
                tot_to_print.append(final_mutation_dict[generation][mutation]['total_cov'])
            except KeyError:  # KeyError should only be raised when mutation not present at generation
                tot_to_print.append("NA")
            try:
                freq_to_print.append(final_mutation_dict[generation][mutation]['frequency'])
            except KeyError:  # KeyError should only be raised when mutation not present at generation
                freq_to_print.append("NA")
        print>>Var, "\t".join(map(str, var_to_print))
        print>>Tot, "\t".join(map(str, tot_to_print))
        print>>Freq, "\t".join(map(str, freq_to_print))
