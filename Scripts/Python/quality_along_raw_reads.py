#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Computes quality scores along raw (unaligned) reads
"""

# --------------------------------------------------------------------------------------------------
# Packages

import sys
import os
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


def get_short_name(long_name):
    L_name = long_name.replace("_", " ").split(" ")
    if len(L_name) == 1:
        return long_name
    L_name[0] = L_name[0][0] + "."
    short_name = " ".join(L_name)
    return short_name


# --------------------------------------------------------------------------------------------------
# Main

if __name__ == "__main__":

    ## Parses arguments
    NUMBER_EXPECTED_ARGUMENTS = 6
    if len(sys.argv) != NUMBER_EXPECTED_ARGUMENTS + 1:
        print(f"ERROR: Wrong number of arguments: {NUMBER_EXPECTED_ARGUMENTS} expected but {len(sys.argv)-1} given.")
        sys.exit(2)
    RAW_READ_DIRNAME, OUTPUT_RAW, OUTPUT_PLOT, FILENAME_COLOR_GC_SPECIES, FILENAME_SPECIES_GROUP, END_LENGTH = sys.argv[1:]
    END_LENGTH = int(END_LENGTH)


    print(f"   - Computes quality scores along raw reads, and at read ends ({END_LENGTH} bases)")

    if os.path.exists(FILENAME_SPECIES_GROUP):
        mode_group = True
    else:
        mode_group = False

    # --- Open output files ---
    # If species are grouped:
    if mode_group:
        # quality along full reads
        output = open(OUTPUT_RAW + "quality_scores_along_reads.txt", "w")
        output.write("\t".join(["Groups"] + list(map(str, np.arange(100)))) + "\n")
        # quality of begining of reads
        output_start = open(OUTPUT_RAW + "quality_scores_along_reads_start.txt", "w")
        output_start.write("\t".join(["Groups"] + list(map(str, np.arange(0, END_LENGTH, 1)))) + "\n")
        # quality of ending of reads
        output_end = open(OUTPUT_RAW + "quality_scores_along_reads_end.txt", "w")
        output_end.write("\t".join(["Groupe"] + list(map(str, np.arange(END_LENGTH, 0, -1)))) + "\n")
    # If species are not grouped:
    else:
        # quality along full reads
        output = open(OUTPUT_RAW + "quality_scores_along_reads.txt", "w")
        output.write("\t".join(["Species"] + list(map(str, np.arange(100)))) + "\n")
        # quality of begining of reads
        output_start = open(OUTPUT_RAW + "quality_scores_along_reads_start.txt", "w")
        output_start.write("\t".join(["Species"] + list(map(str, np.arange(0, END_LENGTH, 1)))) + "\n")
        # quality of ending of reads
        output_end = open(OUTPUT_RAW + "quality_scores_along_reads_end.txt", "w")
        output_end.write("\t".join(["Species"] + list(map(str, np.arange(END_LENGTH, 0, -1)))) + "\n")


    # --- Prepare plot ---
    gs = gridspec.GridSpec(2, 2)
    fig = plt.figure()
    ax = plt.subplot(gs[0, :])
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.95, wspace=0.2, hspace=0.2)
    fig.set_dpi(300.0)
    number_of_split = 100 # split reads into 100 parts (%)


    if mode_group:
        # --- Get color associated to each groups ---
        # and dictionary associating species and group
        dict_group_color = {}
        dict_species_group = {}
        color_file = open(FILENAME_SPECIES_GROUP, "r")
        for line in color_file:
            group_name, sp_name, color = line.rstrip().split("\t")
            sp_name = sp_name.replace(" ", "_")
            sp_name = get_short_name(sp_name)
            if color[0] != "#":
                color = "#" + color
            dict_group_color[group_name] = color
            dict_species_group[sp_name] = group_name
        color_file.close()

    else:
        # --- Get color associated to each species ---
        dict_species_color = {}
        dict_gc_species = {}
        color_file = open(FILENAME_COLOR_GC_SPECIES, "r")
        for line in color_file:
            sp_name, color, gc = line.rstrip().split("\t")
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


    # --- Browse fastq files to get quality scores ---
    if mode_group:
        dict_to_plot = {}
    for raw_read_filename in os.listdir(RAW_READ_DIRNAME):
        if os.path.isdir(raw_read_filename) or ".fastq" not in raw_read_filename:
            continue

        species_name = raw_read_filename.split(".fastq")[0]
        print("  ", species_name)
        species_name_short = get_short_name(species_name)
        if mode_group:
            group_name = dict_species_group[species_name_short]
            if group_name not in dict_to_plot:
                dict_to_plot[group_name] = {}
                for cat in ["full", "start", "end"]:
                    dict_to_plot[group_name][cat] = {}

        # Prepare dictionaries that will store results
        #   Full read: store list of mean quality scores for each % position (0-99)
        dict_position_quality = {}
        for position in range(number_of_split):
            dict_position_quality[position] = [0, 0] # sum of qualities, nb of occurrences
        #   Ends of reads
        dict_position_quality_start, dict_position_quality_end = {}, {}
        for position in range(END_LENGTH):
            dict_position_quality_start[position] = [0, 0]
            dict_position_quality_end[position] = [0, 0]

        # --- Get quality scores for full reads, and "zoom" on reads' ends ---
        raw_read_file = open(RAW_READ_DIRNAME + raw_read_filename, "r")
        while True:
            header = raw_read_file.readline()
            if not header:
                break
            raw_read_file.readline() # skip read line, useless here
            raw_read_file.readline() # skip linking line, useless too
            quality_score_str = raw_read_file.readline().rstrip()

            if len(quality_score_str) < number_of_split: # do not include too short reads
                continue

            # Split quality scores (full read)
            quality_score_list = split(quality_score_str, number_of_split)
            # Store results into dictionary:
            #  for each % position, store mean quality scores for this given position
            for position, quality_scores_part in enumerate(quality_score_list):
                quality_score = get_mean_quality(quality_scores_part)
                dict_position_quality[position][0] += quality_score
                dict_position_quality[position][1] += 1

            # Compute quality on ends of reads
            quality_score_start = quality_score_str[:END_LENGTH]
            quality_score_end = quality_score_str[-END_LENGTH:]
            for position, quality_score in enumerate(quality_score_start):
                dict_position_quality_start[position][0] += ord(quality_score)
                dict_position_quality_start[position][1] += 1
            for position, quality_score in enumerate(quality_score_end):
                dict_position_quality_end[position][0] += ord(quality_score)
                dict_position_quality_end[position][1] += 1

        raw_read_file.close()

        # --- For the given species, compute mean quality scores for each %,
        #     accross each read analysed ---
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


        # If species are grouped, store results for later
        if mode_group:
            for i, position in enumerate(list_positions):
                quality = list_quality_scores[i]
                quality_start = list_quality_scores_start[i]
                quality_end = list_quality_scores_end[i]
                # full read
                if position not in dict_to_plot[group_name]["full"]:
                    dict_to_plot[group_name]["full"][position] = []
                dict_to_plot[group_name]["full"][position] += [quality]
                # read start
                if position not in dict_to_plot[group_name]["start"]:
                    dict_to_plot[group_name]["start"][position] = []
                dict_to_plot[group_name]["start"][position] += [quality_start]
                # read end
                if position not in dict_to_plot[group_name]["end"]:
                    dict_to_plot[group_name]["end"][position] = []
                dict_to_plot[group_name]["end"][position] += [quality_end]


        # Else, plot and write in files now
        else:
            # --- Store results in raw txt file ---
            # for full read
            output.write("\t".join([species_name_short] + list(map(str, list_quality_scores))) + "\n")
            # store results for ends of reads
            output_start.write("\t".join([species_name_short] +
                                         list(map(str, list_quality_scores_start))) + "\n")
            output_end.write("\t".join([species_name_short] +
                                       list(map(str, list_quality_scores_end))) + "\n")


            # Plot results for each species
            plt.plot(list_positions, list_quality_scores, label=species_name_short,
                     color=dict_species_color[species_name_short])
            plt.xticks(ticks=np.arange(0, number_of_split+1, 10),
                       labels=np.arange(0, number_of_split+1, 10))


    if mode_group: # create plot and write results in raw file
        for group_name in dict_to_plot:
            color = dict_group_color[group_name]
            list_positions, list_quality_scores = [], []
            for position in sorted(dict_to_plot[group_name]["full"]):
                list_positions += [position]
                list_quality_scores += [np.mean(dict_to_plot[group_name]["full"][position])]
            plt.plot(list_positions, list_quality_scores, label=group_name, color=color)
            output.write("\t".join([group_name] + list(map(str, list_quality_scores))) + "\n")
        plt.xticks(ticks=np.arange(0, number_of_split+1, 10),
                   labels=np.arange(0, number_of_split+1, 10))
        plt.xlabel("Relative position in read (%)")
        plt.ylabel("Mean quality score")
        plt.legend(ncol=2, title="Groups")
        output.close()


    else: # add labels and legend to the plot
        # Labels and legend
        plt.xlabel("Relative position in read (%)")
        plt.ylabel("Mean quality score")
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ordered_label_list = [elt for elt in list_ordered_species_gc if elt in by_label]
        ordered_label_values = [by_label[k] for k in ordered_label_list]
        plt.legend(ordered_label_values, ordered_label_list, ncol=2, title="Species")

        # Close output files
        output.close()
        output_start.close()
        output_end.close()

    plt.savefig(OUTPUT_PLOT + "quality_scores_along_reads.png")


    # --- Add plots of results for ends of reads ---
    label_step = int(END_LENGTH/80)*10

    if mode_group:
        # Start
        ax = pl.subplot(gs[1, 0])
        list_positions = np.arange(END_LENGTH)
        for group_name in dict_to_plot:
            color = dict_group_color[group_name]
            list_quality_scores = []
            for position in list_positions:
                list_quality_scores += [np.mean(dict_to_plot[group_name]["start"][position])]
            plt.plot(list_positions, list_quality_scores, label=group_name, color=color)
            output_start.write("\t".join([group_name] + list(map(str, list_quality_scores))) + "\n")
        plt.xticks(ticks=np.arange(0, END_LENGTH+1, label_step),
                   labels=np.arange(0, END_LENGTH+1, label_step))
        plt.xlabel("First positions in the read")
        plt.ylabel("Mean quality score")

        # End
        ax = pl.subplot(gs[1, 1])
        for group_name in dict_to_plot:
            color = dict_group_color[group_name]
            list_quality_scores = []
            for position in list_positions:
                list_quality_scores += [np.mean(dict_to_plot[group_name]["end"][position])]
            plt.plot(list_positions, list_quality_scores, label=group_name, color=color)
            output_end.write("\t".join([group_name] + list(map(str, list_quality_scores))) + "\n")
        plt.xticks(ticks=np.arange(0, END_LENGTH+1, label_step),
                   labels=np.arange(END_LENGTH, -1, -label_step))
        plt.xlabel("Last positions in the read")
        plt.ylabel("Mean quality score")

        plt.savefig(OUTPUT_PLOT + "quality_scores_along_reads.png")
        plt.close()
        output_start.close()
        output_end.close()



    else:

        # Start
        ax = pl.subplot(gs[1, 0])
        result_file = open(OUTPUT_RAW + "quality_scores_along_reads_start.txt", "r")
        header = result_file.readline().rstrip()
        list_positions = np.arange(END_LENGTH)
        for line in result_file:
            line = line.rstrip()
            list_qualities_str = line.split("\t")[1:]
            list_qualities = [float(qual) for qual in list_qualities_str]
            species_name = line.split("\t")[0]
            plt.plot(list_positions, list_qualities, color=dict_species_color[species_name])
        result_file.close()
        plt.xticks(ticks=np.arange(0, END_LENGTH+1, label_step),
                   labels=np.arange(0, END_LENGTH+1, label_step))
        plt.xlabel("First positions in the read")
        plt.ylabel("Mean quality score")

        # End
        ax = pl.subplot(gs[1, 1])
        result_file = open(OUTPUT_RAW + "quality_scores_along_reads_end.txt", "r")
        header = result_file.readline().rstrip()
        list_positions = np.arange(END_LENGTH)
        for line in result_file:
            line = line.rstrip()
            list_qualities_str = line.split("\t")[1:]
            list_qualities = [float(qual) for qual in list_qualities_str]
            species_name = line.split("\t")[0]
            plt.plot(list_positions, list_qualities, color=dict_species_color[species_name])
        result_file.close()
        plt.xticks(ticks=np.arange(0, END_LENGTH+1, label_step),
                   labels=np.arange(END_LENGTH, -1, -label_step))
        plt.xlabel("Last positions in the read")
        plt.ylabel("Mean quality score")

        plt.savefig(OUTPUT_PLOT + "quality_scores_along_reads.png")
        plt.close()

