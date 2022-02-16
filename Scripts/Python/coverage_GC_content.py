#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Compute the relative coverage depending on local GC content of reads, for each species
"""

# --------------------------------------------------------------------------------------------------
# Packages

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
PLOT_PARAMS = {'legend.fontsize': 12,
               'legend.title_fontsize': 14,
               'figure.figsize': (14, 9),
               'axes.labelsize': 14,
               'axes.titlesize': 14,
               'xtick.labelsize': 14,
               'ytick.labelsize': 14,
               'legend.borderpad': 0.2,
               'legend.markerscale': 5}

plt.rcParams.update(PLOT_PARAMS)




# --------------------------------------------------------------------------------------------------
# Functions

def get_color(filename):
    """
    Retrieve associated color for each species, from file
    Return a dictionary with species names as keys, and hexadecimal color as values
    Also return a list of ordered species, according to GC content
    """

    # Extract information about GC content and color associated to each species
    dictionary_color = {} # key = species name ; value = color
    dictionary_gc = {} # key = gc content ; value = list of corresponding species name(s)
    file = open(filename, "r")
    for line in file:
        species_name, color, gc = line.rstrip().split("\t")
        gc = float(gc)
        species_name = species_name.replace(" ", "_")
        dictionary_color[species_name] = color
        if gc not in dictionary_gc:
            dictionary_gc[gc] = []
        dictionary_gc[gc] += [species_name]
    file.close()

    # file is supposed to be already sorted in increasing GC content, but sort it here again
    #    to be sure (as user modifications of file are allowed)
    L_ordered_species = [get_short_name(species_name) for gc in sorted(dictionary_gc)
                                                      for species_name in dictionary_gc[gc]]
    return dictionary_color, L_ordered_species



def get_short_name(long_name):
    L_name = long_name.replace("_", " ").split(" ")
    if len(L_name) == 1:
        return long_name
    L_name[0] = L_name[0][0] + "."
    short_name = " ".join(L_name)
    return short_name


def get_gc_content(sequence):
    """
    Returns GC content of a given sequence
    i.e. percentage of (G+C) bases, over all other bases

    """

    nb_C = genome_part.count("C")
    nb_G = genome_part.count("G")
    nb_N = genome_part.count("N")
    if nb_N == WINDOW_SIZE: # ignore, only N bases
        return -1
    return int(round((nb_G + nb_C) / (len(sequence)) * 100))


def get_species_group(filename):
    """
    Returns a dictionary of species and their associated group
    Key = species name
    Value = group name
    Also returns a dictionary of color (value) associated to each group (key)
    """
    dictionary_species = {}
    dictionary_color = {}
    file = open(filename, "r")
    for line in file:
        group_name, species_name, color = line.rstrip().split("\t")
        species_name = get_short_name(species_name.split(".")[0])
        dictionary_species[species_name] = group_name
        dictionary_color[group_name] = color
    file.close()
    return dictionary_species, dictionary_color

def merge_results_groups(dict_unmerged_results, dict_groups):
    """
    From a dictionary of unmerged results (i.e. results for each species separately)
    returns a dictionary of these results, grouped accorded to user choice
    """
    dict_merged_results = {}

    # Gather data
    for species_name in dict_unmerged_results:
        group_name = dict_groups[species_name]
        if group_name not in dict_merged_results:
            dict_merged_results[group_name] = {}
        L_gc, L_coverage = dict_unmerged_results[species_name]
        for i, gc in enumerate(L_gc):
            coverage = L_coverage[i]
            if gc not in dict_merged_results[group_name]:
                dict_merged_results[group_name][gc] = []
            dict_merged_results[group_name][gc] += [coverage]

    # Compute new coverage median
    for group_name in dict_merged_results:
        L_gc, L_coverage = [], []
        for gc in sorted(dict_merged_results[group_name]):
            L_gc += [gc]
            L_coverage += [np.median(dict_merged_results[group_name][gc])]
        dict_merged_results[group_name] = [L_gc, L_coverage]

    return dict_merged_results


# --------------------------------------------------------------------------------------------------
# Main

if __name__ == "__main__":

    ## Parses arguments
    NUMBER_EXPECTED_ARGUMENTS = 9
    if len(sys.argv) != NUMBER_EXPECTED_ARGUMENTS + 1:
        print(f"ERROR: Wrong number of arguments: {NUMBER_EXPECTED_ARGUMENTS} expected but {len(sys.argv)-1} given.")
        sys.exit(2)
    ALN_DIRNAME, GENOME_DIRNAME, FILENAME_COLOR_GC_SPECIES, OUTPUT_RAW, OUTPUT_PLOT, FILE_SPECIES_GROUPS, WINDOW_SIZE, NB_MIN_BASES, NB_MAX_BASES = sys.argv[1:]
    WINDOW_SIZE = int(WINDOW_SIZE)
    NB_MIN_BASES = int(NB_MIN_BASES)
    NB_MAX_BASES = int(NB_MAX_BASES)

    if NB_MAX_BASES == -1:
        NB_MAX_BASES = float("inf")


    # Get information for colors of species depending on their GC content
    dict_species_color, L_ordered_species = get_color(FILENAME_COLOR_GC_SPECIES)

    # Initialize plot and main dictionary that will store data to plot
    dict_species_gc_coverage = {} # key = species name ; value = [gc_content, coverage]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.set_dpi(300.0)
    box = ax.get_position()
    ax.set_position([0.1, 0.1, box.width*0.8, box.height])

    # Parse alignments to get coverage and gc content of reads
    #   and write results in file
    output = open(OUTPUT_RAW + "relative_coverage_gc_content.txt", "w")
    for filename_aln in os.listdir(ALN_DIRNAME):
        if os.path.isdir(ALN_DIRNAME + filename_aln):
            continue

        species_name = filename_aln.split(".txt")[0]
        species_short_name = get_short_name(species_name)

        print("  ", species_name)
        file = open(ALN_DIRNAME + filename_aln, "r")

        # Get coverage for each base in the genome, from alignment file
        #  i.e. get global mean coverage of genome
        dict_coverage_of_genome = {}
        nb_sequenced_bases = 0 # for global expected coverage
        while True:
            header = file.readline()
            if not header:
                break
            genome = file.readline().rstrip()
            read = file.readline().rstrip()
            nb_sequenced_bases += len(read.replace("-", ""))
            pos_in_genome = int(header.split()[6].split("=")[-1]) - 1

            if "- (reverse)" in header:
                step_position = -1
                pos_in_genome += len(genome.replace("-", ""))
            else:
                step_position = 1

            for i in range(len(genome)):
                base_genome = genome[i]
                base_read = read[i]
                if base_read != "-":
                    if pos_in_genome not in dict_coverage_of_genome:
                        dict_coverage_of_genome[pos_in_genome] = 0
                    dict_coverage_of_genome[pos_in_genome] += 1
                if base_genome == "-":
                    continue
                pos_in_genome += step_position
        file.close()

        for i in range(max(dict_coverage_of_genome.keys())):
            if i not in dict_coverage_of_genome:
                dict_coverage_of_genome[i] = 0


        # Browse genome, by sliding window, and compute coverage depending on GC content
        for filename in os.listdir(GENOME_DIRNAME):
            if species_name not in filename:
                continue

            genome = ""
            file = open(GENOME_DIRNAME + filename, "r")
            for line in file:
                if ">" in line:
                    continue
                genome += line.replace("\n", "")
            file.close()

            dict_GC_coverage = {}

            for i in range(len(genome)-WINDOW_SIZE+1):
                genome_part = genome[i: i+WINDOW_SIZE]
                gc_content = get_gc_content(genome_part)
                if gc_content == -1: # means that there are only "N" bases, ignore it
                    continue
                coverage = sum([dict_coverage_of_genome[i] for i in range(i, i+WINDOW_SIZE)]) / WINDOW_SIZE

                if gc_content not in dict_GC_coverage:
                    dict_GC_coverage[gc_content] = []
                dict_GC_coverage[gc_content] += [coverage]

                if i > NB_MAX_BASES:
                    break

            break

        expected_coverage = np.mean([np.mean(elt) for elt in dict_GC_coverage.values()])

        # Plot result
        list_GC = []
        list_mean_coverage = []
        for gc_content in sorted(dict_GC_coverage):
            L_coverage = dict_GC_coverage[gc_content]
            if len(L_coverage) < NB_MIN_BASES:
                continue
            list_GC += [gc_content]
            list_mean_coverage += [np.mean(L_coverage)/expected_coverage]
        ax.plot(list_GC, list_mean_coverage,
                label=species_short_name,
                color=dict_species_color[species_name])

        dict_species_gc_coverage[species_short_name] = [list_GC, list_mean_coverage]

        # Write
        output.write(species_name + "(GC then coverage)\n")
        output.write("\t".join(map(str, list_GC)) + "\n")
        output.write("\t".join(map(str, list_mean_coverage)) + "\n\n")

    output.close()

    # Plot theoretical line of unbiased coverage depending on GC content
    unbiased_x = list(np.arange(0, 101))
    unbiased_y = [1] * len(unbiased_x)
    ax.plot(unbiased_x, unbiased_y, "--", color="k")


    # Reorder legend labels
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ordered_label_list = [elt for elt in L_ordered_species if elt in by_label]
    ordered_label_values = [by_label[k] for k in ordered_label_list]
    #plt.legend(ordered_label_values, ordered_label_list, ncol=3, title="Species:")
    ax.legend(ordered_label_values, ordered_label_list, title="Species", bbox_to_anchor=(1.35,1))


    plt.xlim(-1, 101)
    plt.xticks(np.arange(0, 101, 10))

    plt.xlabel(f"% GC ({WINDOW_SIZE}-base window)")
    plt.ylabel("Relative coverage")

    plt.savefig(OUTPUT_PLOT + "relative_coverage_gc_content.png")
    plt.close()




    # --- Plot another graph: same but with data grouped ---
    #  (only if groups were defined by user)
    if not os.path.exists(FILE_SPECIES_GROUPS):
        sys.exit(0)
    dict_species_group, dict_color_group = get_species_group(FILE_SPECIES_GROUPS)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.set_dpi(300.0)
    box = ax.get_position()
    ax.set_position([0.1, 0.1, box.width*0.8, box.height])
    dict_species_gc_coverage_grouped = merge_results_groups(dict_species_gc_coverage,
                                                            dict_species_group)
    for group_name in dict_species_gc_coverage_grouped:
        list_GC, list_mean_coverage = dict_species_gc_coverage_grouped[group_name]
        plt.plot(list_GC, list_mean_coverage, label=group_name,
                color=dict_color_group[group_name])
    # theoretical unbiased coverage
    unbiased_x = list(np.arange(0, 101))
    unbiased_y = [1] * len(unbiased_x)
    ax.plot(unbiased_x, unbiased_y, "--", color="k")
    ax.legend(title="Groups", bbox_to_anchor=(1.35,1))
    # axis labels
    plt.xlim(-1, 101)
    plt.xticks(np.arange(0, 101, 10))
    plt.xlabel(f"% GC ({WINDOW_SIZE}-base window)")
    plt.ylabel("Relative coverage")
    # save
    plt.savefig(OUTPUT_PLOT + "relative_coverage_gc_content_grouped.png")
    plt.close()



