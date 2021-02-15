#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Computes quality scores along raw (unaligned) reads
"""

# --------------------------------------------------------------------------------------------------
# Packages

import sys
import os
import ast
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pylab as pl
import matplotlib.gridspec as gridspec

PARAMS = {'legend.fontsize': 14,
          'legend.title_fontsize': 14,
          'legend.labelspacing': 0.1,
          'legend.borderpad': 0.3,
          'legend.handletextpad': 0.1,
          'figure.figsize': (14, 9),
          'axes.labelsize': 18,
          'axes.titlesize': 20,
          'xtick.labelsize': 18,
          'ytick.labelsize': 18}
plt.rcParams.update(PARAMS)


# --------------------------------------------------------------------------------------------------
# Functions


def split(string: str, n: int):
    """ Split string into n (almost) even elements

    Parameters
    ----------
    string: the string to split
    n:      desired number of splits

    Returns
    -------
    L: a list of length n, each element containing (almost) same number of characters
    """
    quotien, rest = divmod(len(string), n)
    L = [string[i * quotien + min(i, rest):(i + 1) * quotien + min(i + 1, rest)] for i in range(n)]
    return L


def get_mean_quality(quality_str):
    """ Return mean quality of given string of qualities (-33 to fit Nanopore quality scores)
    """
    quality_score = 0
    for letter in quality_str:
        quality_score += ord(letter)
    mean_quality_score = quality_score / len(quality_str)
    mean_quality_score -= 33
    return mean_quality_score

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

    ## Parses arguments

    if len(sys.argv) != 6:
        print(f"ERROR: Wrong number of arguments: 5 expected but {len(sys.argv)-1} given.")
        sys.exit(2)
    RAW_READ_DIRNAME = sys.argv[1]
    ALN_DIRNAME = sys.argv[2]
    OUTPUT_RAW = sys.argv[3]
    OUTPUT_PLOT = sys.argv[4]
    FILENAME_COLOR_GC_SPECIES = sys.argv[5]

    if RAW_READ_DIRNAME[-1] != "/":
        RAW_READ_DIRNAME += "/"
    if ALN_DIRNAME[-1] != "/":
        ALN_DIRNAME += "/"
    if OUTPUT_PLOT[-1] != "/":
        OUTPUT_PLOT += "/"
    if OUTPUT_RAW[-1] != "/":
        OUTPUT_RAW += "/"

    output = open(OUTPUT_RAW + "quality_scores_along_reads.txt", "w")
    output.write("\t".join(["Species"] + list(map(str, np.arange(100)))) + "\n")
    output_start = open(OUTPUT_RAW + "quality_scores_along_reads_start.txt", "w")
    output_start.write("\t".join(["Species"] + list(map(str, np.arange(0, 200, 1)))) + "\n")
    output_end = open(OUTPUT_RAW + "quality_scores_along_reads_end.txt", "w")
    output_end.write("\t".join(["Species"] + list(map(str, np.arange(200, 0, -1)))) + "\n")

    number_of_split = 100

    gs = gridspec.GridSpec(2, 2)

    fig = plt.figure()
    ax = plt.subplot(gs[0, :])
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.95, wspace=0.2, hspace=0.2)
    fig.set_dpi(300.0)



    for raw_read_filename in os.listdir(RAW_READ_DIRNAME):
        if ".fastq" not in raw_read_filename:
            continue

        if os.path.isdir(raw_read_filename) or ".fastq" not in raw_read_filename:
            continue

        print(raw_read_filename)

        species_name = raw_read_filename.split(".fastq")[0]
        species_name_short = get_short_name(species_name)
        dict_position_quality = {} # Store list of mean quality scores for each % position (0-99)
        for position in range(number_of_split):
            dict_position_quality[position] = [0, 0] # sum of qualities, nb of occurrences

        # Zoom on both ends of reads
        dict_position_quality_start, dict_position_quality_end = {}, {}
        for position in range(200):
            dict_position_quality_start[position] = [0, 0]
            dict_position_quality_end[position] = [0, 0]

        raw_read_file = open(RAW_READ_DIRNAME + raw_read_filename, "r")
        while True:
            header = raw_read_file.readline()
            if not header:
                break
            raw_read_file.readline() # skip read line, useless here
            raw_read_file.readline() # skip useless linking line
            quality_score_str = raw_read_file.readline().rstrip()

            if len(quality_score_str) < number_of_split:
                continue

            # Split quality scores into 100 parts
            quality_score_list = split(quality_score_str, number_of_split)
            # Store results into dictionary:
            #  for each % position, store mean quality scores for this given position
            for position, quality_scores_part in enumerate(quality_score_list):
                quality_score = get_mean_quality(quality_scores_part)
                dict_position_quality[position][0] += quality_score
                dict_position_quality[position][1] += 1


            # Zoom on both ends of reads
            quality_score_start = quality_score_str[:200]
            quality_score_end = quality_score_str[-200:]
            for position, quality_score in enumerate(quality_score_start):
                dict_position_quality_start[position][0] += ord(quality_score)
                dict_position_quality_start[position][1] += 1
            for position, quality_score in enumerate(quality_score_end):
                dict_position_quality_end[position][0] += ord(quality_score)
                dict_position_quality_end[position][1] += 1

        raw_read_file.close()

        # For the given species, compute mean quality scores for each %, accross each read analysed
        list_positions, list_quality_scores = [], []
        for position in sorted(dict_position_quality):
            list_positions += [position]
            sum_qualities, nb_occurrences = dict_position_quality[position]
            list_quality_scores += [sum_qualities / nb_occurrences]

        list_quality_scores_start, list_quality_scores_end = [], []
        for position in sorted(dict_position_quality_start):
            sum_qualities, nb_occurrences = dict_position_quality_start[position]
            list_quality_scores_start += [sum_qualities / nb_occurrences]
            sum_qualities, nb_occurrences = dict_position_quality_end[position]
            list_quality_scores_end += [sum_qualities / nb_occurrences]

        #  reorder legend's labels
        dict_species_color = {}
        dict_gc_species = {}
        color_file = open(FILENAME_COLOR_GC_SPECIES, "r")
        for line in color_file:
            sp_name, color, gc = line.rstrip().split(" ; ")
            sp_name = sp_name.replace(" ", "_")
            sp_name = get_short_name(sp_name)
            gc = float(gc)
            if color[0] != "#":
                color = "#" + color
            dict_species_color[sp_name] = color
            if gc not in dict_gc_species:
                dict_gc_species[gc] = []
            dict_gc_species[gc] += [sp_name]
        color_file.close()
        list_ordered_species_gc = []
        for gc in sorted(dict_gc_species):
            for elt in dict_gc_species[gc]:
                list_ordered_species_gc += [elt]

        #      plot results for each species
        plt.plot(list_positions, list_quality_scores, label=species_name_short,
                 color=dict_species_color[species_name_short])
        plt.xticks(ticks=np.arange(0, number_of_split+1, 10),
                   labels=np.arange(0, number_of_split+1, 10))

        plt.xlabel("Relative position in read (%)")
        plt.ylabel("Mean quality score")


        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ordered_label_list = [elt for elt in list_ordered_species_gc if elt in by_label]
        ordered_label_values = [by_label[k] for k in ordered_label_list]

        plt.legend(ordered_label_values, ordered_label_list, ncol=3, title="Species")

        plt.savefig(OUTPUT_PLOT + "quality_scores_along_reads.png")

        # Store results in raw txt file
        output.write("\t".join([species_name_short] + list(map(str, list_quality_scores))) + "\n")
        # store results for ends of reads
        output_start.write("\t".join([species_name_short] +
                                     list(map(str, list_quality_scores_start))) + "\n")
        output_end.write("\t".join([species_name_short] +
                                   list(map(str, list_quality_scores_end))) + "\n")




    # Close files
    output.close()
    output_start.close()
    output_end.close()


    # Add plots of results for ends of reads

    #  reorder legend's labels

    plt.legend(ordered_label_values, ordered_label_list, ncol=3, title="Species")


    # start
    ax = pl.subplot(gs[1, 0])
    pl.plot([0,1])

    result_file = open(OUTPUT_RAW + "quality_scores_along_reads_start.txt", "r")

    header = result_file.readline().rstrip()
    list_positions = np.arange(200)
    for line in result_file:
        line = line.rstrip()
        list_qualities_str = line.split("\t")[1:]
        list_qualities = [float(qual) for qual in list_qualities_str]
        species_name = line.split("\t")[0]
        plt.plot(list_positions, list_qualities, color=dict_species_color[species_name])
    result_file.close()
    plt.xticks(ticks=np.arange(0, 201, 20), labels=np.arange(0, 201, 20))
    plt.ylim(33, 65)
    plt.xlabel("n first positions in the read")
    plt.ylabel("Mean quality score")

    # end
    ax = pl.subplot(gs[1, 1])
    pl.plot([0,1])
    result_file = open(OUTPUT_RAW + "quality_scores_along_reads_end.txt", "r")
    header = result_file.readline().rstrip()
    list_positions = np.arange(200)
    for line in result_file:
        line = line.rstrip()
        list_qualities_str = line.split("\t")[1:]
        list_qualities = [float(qual) for qual in list_qualities_str]
        species_name = line.split("\t")[0]
        plt.plot(list_positions, list_qualities, color=dict_species_color[species_name])
    result_file.close()
    plt.xticks(ticks=np.arange(0, 201, 20), labels=np.arange(200, -1, -20))
    plt.ylim(33, 65)
    plt.xlabel("n last positions in the read")
    plt.ylabel("Mean quality score")

    plt.savefig(OUTPUT_PLOT + "quality_scores_along_reads.png")
    plt.close()
