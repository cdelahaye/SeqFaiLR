#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Compute the relative coverage depending on local GC content of reads, for each bacterial species
"""

# --------------------------------------------------------------------------------------------------
# Packages

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
PLOT_PARAMS = {'legend.fontsize': 15,
               'legend.title_fontsize': 17,
               'legend.labelspacing':0.2,
               'figure.figsize': (14, 9),
               'axes.labelsize': 18,
               'axes.titlesize': 18,
               'xtick.labelsize': 18,
               'ytick.labelsize': 18,
               'legend.borderpad': 0.2,
               'legend.markerscale': 5}

plt.rcParams.update(PLOT_PARAMS)




# --------------------------------------------------------------------------------------------------
# Functions

def get_short_name(long_species_name: str):
    """ Returns a short version of species name
         and replace underscores with spaces
    For example: Staphylococcus_thermophilus_CNRZ1066 -> S. thermophilus CNRZ1066
    """
    list_long_species_name = long_species_name.split("_")
    short_species_name = list_long_species_name[0][0] + ". "
    short_species_name += " ".join(list_long_species_name[1:])
    return short_species_name



# --------------------------------------------------------------------------------------------------
# Main

if __name__ == "__main__":

    # Parse arguments
    if len(sys.argv) != 5:
        print(f"ERROR: Wrong number of arguments: 4 expected but {len(sys.argv)-1} given.")
        sys.exit(2)
    ALN_DIRNAME = sys.argv[1]
    GENOME_DIRNAME = sys.argv[2]
    FILENAME_SPECIES_GC_COLOR = sys.argv[3]
    OUTPUT_PLOT = sys.argv[4]

    if ALN_DIRNAME[-1] != "/":
        ALN_DIRNAME += "/"
    if GENOME_DIRNAME[-1] != "/":
        GENOME_DIRNAME += "/"
    if OUTPUT_PLOT[-1] != "/":
        OUTPUT_PLOT += "/"


    WINDOW_SIZE = 100

    # Get information for colors of species depending on their GC content
    dict_species_color = {}
    dict_gc_species = {}
    color_file = open(FILENAME_SPECIES_GC_COLOR, "r")
    for line in color_file:
        species_name, color, gc = line.rstrip().split(" ; ")
        species_name = species_name.replace(" ", "_")
        gc = float(gc)
        if color[0] != "#":
            color = "#" + color
        dict_species_color[species_name] = color
        if gc not in dict_gc_species:
            dict_gc_species[gc] = []
        dict_gc_species[gc] += [species_name]
    color_file.close()
    list_ordered_species_gc = []
    for gc in sorted(dict_gc_species):
        for species in dict_gc_species[gc]:
            list_ordered_species_gc += [get_short_name(species)]


    dict_store_results = {}
    fig, ax = plt.subplots()

    for filename_aln in os.listdir(ALN_DIRNAME):

        if os.path.isdir(ALN_DIRNAME + filename_aln):
            continue

        species_name = filename_aln.split(".txt")[0]
        species_short_name = get_short_name(species_name)

        print(species_name)

        file = open(ALN_DIRNAME + filename_aln, "r")
        dict_coverage_of_genome = {}


        nb_sequenced_bases = 0 # for global expected coverage


        # Get coverage for each base in the genome, from alignment file
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

        expected_coverage = nb_sequenced_bases / max(dict_coverage_of_genome.keys())

        # Browse genome, by sliding window, and compute
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

                nb_A = genome_part.count("A")
                nb_C = genome_part.count("C")
                nb_G = genome_part.count("G")
                nb_T = genome_part.count("T")
                gc_content = int(round((nb_G + nb_C) / (nb_A + nb_C + nb_G + nb_T) * 100))

                coverage = sum([dict_coverage_of_genome[i] for i in range(i, i+WINDOW_SIZE)]) / WINDOW_SIZE

                if gc_content not in dict_GC_coverage:
                    dict_GC_coverage[gc_content] = []
                dict_GC_coverage[gc_content] += [coverage]

            break

        expected_coverage = np.mean([np.mean(elt) for elt in dict_GC_coverage.values()])

        # Plot result
        list_GC = []
        list_mean_coverage = []
        for gc_content in sorted(dict_GC_coverage):
            L_coverage = dict_GC_coverage[gc_content]
            if len(L_coverage) < 1000:
                continue
            list_GC += [gc_content]
            list_mean_coverage += [np.mean(L_coverage)/expected_coverage]
        ax.plot(list_GC, list_mean_coverage,
                label=species_short_name,
                color=dict_species_color[species_name])

        dict_store_results[species_short_name] = [list_GC, list_mean_coverage]


    unbiased_x = list(np.arange(0, 101))
    unbiased_y = [1] * len(unbiased_x)
    ax.plot(unbiased_x, unbiased_y, "--", color="k")


    # Reorder legend labels
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ordered_label_list = [elt for elt in list_ordered_species_gc if elt in by_label]
    ordered_label_values = [by_label[k] for k in ordered_label_list]
    plt.legend(ordered_label_values, ordered_label_list, ncol=3, title="Species:")



    plt.xlim(0, 100)
    plt.xticks(np.arange(0, 100, 10))
    plt.ylim(0.7, 1.3)
    plt.yticks(np.arange(0.7, 1.4, 0.1))

    plt.xlabel("% GC (100-base window)")
    plt.ylabel("Relative coverage")

    plt.savefig(OUTPUT_PLOT + "relative_coverage_gc_content.png")
