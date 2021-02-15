#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
For each category of gc content (low / high gc content for bacteria + human), computes
mean error rates (total and detailled) depending on the GC content of the read
"""


# --------------------------------------------------------------------------------------------------
# Packages

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

params = {'legend.fontsize': 16,
          'legend.title_fontsize': 16,
          'legend.labelspacing': 0.1,
          'legend.borderpad': 0.3,
          'legend.columnspacing': 1,
          'legend.handletextpad': 0.5,
          'legend.handlelength': 1.2,
          'figure.figsize': (12, 8),
          'axes.labelsize': 16,
          'axes.titlesize': 16,
          'xtick.labelsize': 16,
          'ytick.labelsize': 16}
plt.rcParams.update(params)

# --------------------------------------------------------------------------------------------------
# Functions

def get_error_rates() :
    """ Computes error rates (total, mismatch, insertion and deletion) for a given alignment
    between a genome and and read
    """
    mismatches, insertions, deletions = 0, 0, 0
    length = len(genome.replace("N", ""))
    for i, base_genome in enumerate(genome) :
        base_read = read[i]
        if base_genome == "N" :
            continue
        elif base_genome == "-" :
            insertions += 1
        elif base_read == "-" :
            deletions += 1
        elif base_genome != base_read :
            mismatches += 1
    total = round((mismatches + insertions + deletions) / length * 100, 1)
    mismatches = round(mismatches / length * 100, 2)
    insertions = round(insertions / length * 100, 2)
    deletions = round(deletions / length * 100, 2)
    return [total, mismatches, insertions, deletions]


def get_gc_content_of_read():
    """ Computes the GC content of a given read
    """
    r = read.replace("-", "")
    gc_content = round((r.count("G") + r.count("C")) / len(r) * 100)
    return gc_content


# --------------------------------------------------------------------------------------------------
# Main

if __name__ == "__main__":


    ## Parses arguments
    if len(sys.argv) != 5:
        print(f"ERROR: Wrong number of arguments: 4 expected but {len(sys.argv)-1} given.")
        sys.exit(2)
    ALN_EXPL_DIRNAME = sys.argv[1]
    OUTPUT_RAW = sys.argv[2]
    OUTPUT_PLOT = sys.argv[3]
    FILENAME_SPECIES_GC = sys.argv[4]


    if ALN_EXPL_DIRNAME[-1] != "/":
        ALN_EXPL_DIRNAME += "/"
    if OUTPUT_PLOT[-1] != "/":
        OUTPUT_PLOT += "/"
    if OUTPUT_RAW[-1] != "/":
        OUTPUT_RAW += "/"


    # Get GC category for each species
    dict_species_gc_category = {}
    file = open(FILENAME_SPECIES_GC, "r")
    for line in file:
        species_name, _, gc = line.rstrip().split(" ; ")
        species_name = species_name.replace(" ", "_")
        gc = float(gc)
        if "Human" in species_name:
            category = "human"
        elif gc < 50:
            category = "low GC"
        else:
            category = "high GC"
        dict_species_gc_category[species_name] = category
    file.close()


    # Initiate dictionary that will store, for each GC category, the error rates depending on
    #   GC content of reads
    dict_error_rates = {}
    list_error_type = ["total", "mismatch", "insertion", "deletion"]


    # Parse alignments to get error rates
    for aln_filename in os.listdir(ALN_EXPL_DIRNAME):

        species_name = aln_filename.split(".txt")[0]
        gc_category = dict_species_gc_category[species_name]
        print(species_name, gc_category)
        file = open(ALN_EXPL_DIRNAME + aln_filename, "r")

        dict_error_rates[species_name] = {}
        for error_type in list_error_type:
            dict_error_rates[species_name][error_type] = {}

        aln_count = 0 # count number of alignments processed
        while True:
            header = file.readline()
            if not header:
                break
            genome = file.readline()
            read = file.readline()
            aln_count += 1
            #if aln_count > 500000: # bound number of analysed alignments
            #    break
            if aln_count % 1000 == 0:
                print(f"{aln_count} alignment analysed")
            gc_content_read = get_gc_content_of_read()
            list_error_rates = get_error_rates()

            for i, error_rate in enumerate(list_error_rates):
                error_type = list_error_type[i]
                if gc_content_read not in dict_error_rates[species_name][error_type]:
                    dict_error_rates[species_name][error_type][gc_content_read] = [0, 0]
                dict_error_rates[species_name][error_type][gc_content_read][0] += error_rate
                dict_error_rates[species_name][error_type][gc_content_read][1] += 1
        file.close()

        # plot for bacteria
        dict_gc_read_errors = {}
        threshold = 100 # minimum occurrences to plot
        for species_name in dict_error_rates:
            gc_category = dict_species_gc_category[species_name]
            if gc_category not in ["low GC", "high GC"]:
                continue
            for gc_content_read in sorted(dict_error_rates[species_name]["mismatch"]):
                total_mismatches, occ_mismatches = dict_error_rates[species_name]["mismatch"][gc_content_read]
                total_insertions, occ_insertions = dict_error_rates[species_name]["insertion"][gc_content_read]
                total_deletions, occ_deletions = dict_error_rates[species_name]["deletion"][gc_content_read]
                if gc_content_read not in dict_gc_read_errors:
                    dict_gc_read_errors[gc_content_read] = [[], [], [], 0]
                dict_gc_read_errors[gc_content_read][0] += [total_mismatches / occ_mismatches]
                dict_gc_read_errors[gc_content_read][1] += [total_insertions / occ_insertions]
                dict_gc_read_errors[gc_content_read][2] += [total_deletions / occ_deletions]
                dict_gc_read_errors[gc_content_read][3] += occ_deletions


        list_mismatches, list_insertions, list_deletions = [], [], []
        list_gc_read = []
        for gc_content in sorted(dict_gc_read_errors):
            if dict_gc_read_errors[gc_content][3] < threshold:
                continue
            list_gc_read += [gc_content]
            list_mismatches += [np.median(dict_gc_read_errors[gc_content][0])]
            list_insertions += [np.median(dict_gc_read_errors[gc_content][1])]
            list_deletions += [np.median(dict_gc_read_errors[gc_content][2])]

        plt.plot(list_gc_read, list_mismatches, label="Mismatch", color="#DDAA33")
        plt.plot(list_gc_read, list_insertions, label="Insertion", color="#BB5566")
        plt.plot(list_gc_read, list_deletions, label="Deletion", color="#004488")
        #  total for high and low gc bacteria
        dict_gc_read_errors_low, dict_gc_read_errors_high = {}, {}
        for species_name in dict_error_rates:
            gc_category = dict_species_gc_category[species_name]
            for gc_content_read in sorted(dict_error_rates[species_name]["total"]):
                total_err, occ_total_err = dict_error_rates[species_name]["total"][gc_content_read]
                if occ_total_err < threshold:
                    continue
                if gc_category == "low GC":
                    if gc_content_read not in dict_gc_read_errors_low:
                        dict_gc_read_errors_low[gc_content_read] = [[], 0]
                    dict_gc_read_errors_low[gc_content_read][0] += [total_err / occ_total_err]
                    dict_gc_read_errors_low[gc_content_read][1] += occ_total_err
                elif gc_category == "high GC":
                    if gc_content_read not in dict_gc_read_errors_high:
                        dict_gc_read_errors_high[gc_content_read] = [[], 0]
                    dict_gc_read_errors_high[gc_content_read][0] += [total_err / occ_total_err]
                    dict_gc_read_errors_high[gc_content_read][1] += occ_total_err
        list_gc_low, list_gc_high = [], []
        list_total_err_low, list_total_err_high = [], []
        for gc_content_read in sorted(dict_gc_read_errors_low):
            if dict_gc_read_errors_low[gc_content_read][1] < threshold:
                continue
            list_gc_low += [gc_content_read]
            list_total_err_low += [np.median(dict_gc_read_errors_low[gc_content_read][0])]
        for gc_content_read in sorted(dict_gc_read_errors_high):
            if dict_gc_read_errors_high[gc_content_read][1] < threshold:
                continue
            list_gc_high += [gc_content_read]
            list_total_err_high += [np.median(dict_gc_read_errors_high[gc_content_read][0])]
        plt.plot(list_gc_low, list_total_err_low, "-.", linewidth=2.5, color="black")
        plt.plot(list_gc_high, list_total_err_high, "--", linewidth=2.5, color="black")
        plt.ylim(0, 10)
        plt.xlim(25, 70)
        plt.xticks(np.arange(25, 71, 5))

        plt.xlabel('GC content of sequenced read (%)')
        plt.ylabel('Error rates (%)')
        # Legend
        deletion_label = Line2D([0], [0], color="#004488", linewidth=2, linestyle='-')
        mismatch_label = Line2D([0], [0], color="#DDAA33", linewidth=2, linestyle='-')
        insertion_label = Line2D([0], [0], color="#BB5566", linewidth=2, linestyle='-')
        tot_err_low_label = Line2D([0], [0], color="k", linewidth=2, linestyle='-.')
        tot_err_high_label = Line2D([0], [0], color="k", linewidth=2, linestyle='--')
        L_lines = [deletion_label, mismatch_label, insertion_label,
                   tot_err_low_label, tot_err_high_label]
        labels = ["Deletions", "Mismatch", "Insertion", "Total (low GC)", "Total (high GC)"]
        plt.legend(L_lines, labels, loc=2, ncol=2)
        plt.savefig(OUTPUT_PLOT + "error_rate_gc_read_bacteria.png")
        plt.close()

        # Write bacterial results in raw file
        file = open(OUTPUT_RAW + "error_rate_gc_read_bacteria.txt", "w")
        file.write("---\n")
        file.write("GC content of reads and total error rates for low GC bacteria\n")
        file.write(f"{list_gc_low}\n")
        file.write(f"{list_total_err_low}\n")
        file.write("---\n")
        file.write("GC content of reads and total error rates for high GC bacteria\n")
        file.write(f"{list_gc_high}\n")
        file.write(f"{list_total_err_high}\n")
        file.write("---\n")
        file.write("GC content of reads and mismatch/insertion/deletion error rates for bacteria\n")
        file.write(f"{list_gc_read}\n")
        file.write(f"{list_mismatches}\n{list_insertions}\n{list_deletions}\n")
        file.close()


        # plot for human
        dict_gc_read_errors = {}
        threshold = 1000
        for species_name in dict_error_rates:
            gc_category = dict_species_gc_category[species_name]
            if gc_category not in ["human"]:
                continue
            for gc_content_read in sorted(dict_error_rates[species_name]["mismatch"]):
                total_mismatches, occ_mismatches = dict_error_rates[species_name]["mismatch"][gc_content_read]
                total_insertions, occ_insertions = dict_error_rates[species_name]["insertion"][gc_content_read]
                total_deletions, occ_deletions = dict_error_rates[species_name]["deletion"][gc_content_read]
                if occ_mismatches < threshold:
                    continue
                if gc_content_read not in dict_gc_read_errors:
                    dict_gc_read_errors[gc_content_read] = [[], [], [], 0]
                dict_gc_read_errors[gc_content_read][0] += [total_mismatches / occ_mismatches]
                dict_gc_read_errors[gc_content_read][1] += [total_insertions / occ_insertions]
                dict_gc_read_errors[gc_content_read][2] += [total_deletions / occ_deletions]
                dict_gc_read_errors[gc_content_read][3] += occ_deletions

        list_mismatches, list_insertions, list_deletions = [], [], []
        list_gc_read = []
        for gc_content in sorted(dict_gc_read_errors):
            if dict_gc_read_errors[gc_content][3] < threshold:
                continue
            list_gc_read += [gc_content]
            list_mismatches += [np.median(dict_gc_read_errors[gc_content][0])]
            list_insertions += [np.median(dict_gc_read_errors[gc_content][1])]
            list_deletions += [np.median(dict_gc_read_errors[gc_content][2])]

        plt.plot(list_gc_read, list_mismatches, label="Mismatch", color="#DDAA33")
        plt.plot(list_gc_read, list_insertions, label="Insertion", color="#BB5566")
        plt.plot(list_gc_read, list_deletions, label="Deletion", color="#004488")
        #  total
        dict_gc_read_errors_human = {}
        for species_name in dict_error_rates:
            gc_category = dict_species_gc_category[species_name]
            if gc_category != "human":
                continue
            for gc_content_read in sorted(dict_error_rates[species_name]["total"]):
                total_err, occ_total_err = dict_error_rates[species_name]["total"][gc_content_read]
                if occ_total_err < threshold:
                    continue
                if gc_content_read not in dict_gc_read_errors_human:
                    dict_gc_read_errors_human[gc_content_read] = [[], 0]
                dict_gc_read_errors_human[gc_content_read][0] += [total_err / occ_total_err]
                dict_gc_read_errors_human[gc_content_read][1] += occ_total_err
        list_gc_read, list_total_err_human = [], []
        for gc_content_read in sorted(dict_gc_read_errors_human):
            list_gc_read += [gc_content_read]
            list_total_err_human += [np.median(dict_gc_read_errors_human[gc_content_read][0])]
        plt.plot(list_gc_read, list_total_err_human, "--", linewidth=2.5, color="black")

        plt.ylim(0, 10)
        plt.xlim(25, 70)
        plt.xticks(np.arange(25, 71, 5))
        plt.xlabel('GC content of sequenced read (%)')
        plt.ylabel('Error rates (%)')

        # Legend
        deletion_label = Line2D([0], [0], color="#004488", linewidth=2, linestyle='-')
        mismatch_label = Line2D([0], [0], color="#DDAA33", linewidth=2, linestyle='-')
        insertion_label = Line2D([0], [0], color="#BB5566", linewidth=2, linestyle='-')
        tot_err_label = Line2D([0], [0], color="k", linewidth=2, linestyle='--')

        L_lines = [deletion_label, mismatch_label, insertion_label, tot_err_label]
        labels = ["Deletions", "Mismatch", "Insertion", "Total"]
        plt.legend(L_lines, labels, loc=2)

        plt.savefig(OUTPUT_PLOT + "error_rate_gc_read_human.png")
        plt.close()


        # Write human results in raw file
        file = open(OUTPUT_RAW + "error_rate_gc_read_human.txt", "w")
        file.write("---\n")

        file.write("---\n")
        file.write("GC content of reads and total/mismatch/insertion/deletion error rates for human\n")
        file.write(f"{list_gc_read}\n")
        file.write(f"{list_total_err_human}\n{list_mismatches}\n{list_insertions}\n{list_deletions}\n")
        file.close()
