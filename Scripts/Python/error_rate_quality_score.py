#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Computes error rates depending on quality score, at read- and base-level, compared to the
expected distribution of Phred scores
"""

# --------------------------------------------------------------------------------------------------
# Packages

import os
import sys
import ast
import matplotlib.pyplot as plt
import numpy as np

params = {'legend.fontsize': 12,
          'legend.title_fontsize': 14,
          'figure.figsize': (14, 9),
          'axes.labelsize': 14,
          'axes.titlesize': 14,
          'xtick.labelsize': 14,
          'ytick.labelsize': 14}
plt.rcParams.update(params)


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


def get_error_rate(genome, read):
    """ Computes error rate between a given alignment of a read against a genome
    """
    error_count = 0
    for i, base_genome in enumerate(genome):
        base_read = read[i]
        if base_read != base_genome:
            error_count += 1
    error_rate = round(error_count / len(genome) * 100, 4)
    return error_rate


def get_quality_read(species_name, read_file, read_id, soft_clips):
    """ Computes mean quality of the current read
    """
    while True:
        header = read_file.readline()
        if not header:
            read_file.close()
            read_file = open(RAW_READ_DIRNAME + species_name + ".fastq", "r")
            continue
        read_file.readline() # skip read sequence
        read_file.readline() # skip useless line
        quality_str = read_file.readline().rstrip()
        if header.split()[0][1:] == read_id:
            break
    quality_score = 0
    soft_clip_start, soft_clip_end = soft_clips
    quality_str = quality_str[soft_clip_start:]
    if soft_clip_end > 0:
        quality_str = quality_str[: -soft_clip_end]
    for char in quality_str:
        quality_score += ord(char)
    quality_score_mean = round(quality_score / len(quality_str), 1) - 33
    return quality_score_mean, read_file


def compute_err_rate_quality_window(dictionary, species_name, genome, read, read_file, read_id, soft_clips):
    """ Computes error rates and gc content for current alignment,
    for WINDOW_LENGTH-bases sliding windows along read
    Updates dictionary storing results
    """

    # Retrieve quality from fastq file
    while True:
        header = read_file.readline()
        if not header:
            read_file.close()
            read_file = open(RAW_READ_DIRNAME + species_name + ".fastq", "r")
            continue
        read_file.readline() # skip read sequence
        read_file.readline() # skip useless line
        quality_str = read_file.readline().rstrip()
        if header.split()[0][1:] == read_id:
            break

    # Compute error rates and quality score for sliding windows
    start_pos, end_pos = 0, WINDOW_LENGTH
    # Initiate part of the quality string
    soft_clip_start, soft_clip_end = soft_clips
    quality_str = quality_str[soft_clip_start:]
    if soft_clip_end > 0:
        quality_str = quality_str[: -soft_clip_end]
    quality_str_part = quality_str[:WINDOW_LENGTH]
    quality_score = 0
    for char in quality_str_part:
        quality_score += ord(char) - 33
    while True:

        # Initiate read and genome portions
        r = read[start_pos: end_pos]
        while len(r.replace("-", "")) < WINDOW_LENGTH and end_pos < len(read)-1:
            r += read[end_pos]
            end_pos += 1
            if end_pos >= len(read)-1:
                return read_file
        if len(r.replace("-", "")) < WINDOW_LENGTH:
                return read_file
        g = genome[start_pos: end_pos]

        # Compute mean error rate
        error_count = 0
        for i, base_g in enumerate(g):
            base_r = r[i]
            if base_g != base_r:
                error_count += 1
        error_rate = round(error_count / len(g) * 100, 2)

        # Compute mean quality score
        quality_score_mean = round(quality_score / WINDOW_LENGTH, 1)

        # Update dictionary storing results
        if quality_score_mean not in dictionary:
            dictionary[quality_score_mean] = [0, 0]
        dictionary[quality_score_mean][0] += error_rate
        dictionary[quality_score_mean][1] += 1

        # Update quality portion and score
        quality_score = quality_score - (ord(quality_str[0]) - 33)
        quality_str = quality_str[1:]
        quality_str_part = quality_str[0: WINDOW_LENGTH]
        quality_score = quality_score + (ord(quality_str_part[-1]) - 33)

        # Update start positions
        start_pos += 1
        while read[start_pos] == "-":
            start_pos += 1


def plot_quality_error_for_reads():
    """
    Computes error rate depending on quality scores, at read level
    And plot results
    In case species were not grouped
    """

    # File to save raw results
    output_raw = open(OUTPUT_RAW + "error_rate_quality_score_reads.txt", "w")
    output_raw.write("Species_name\tMean quality score reads\tMean error rate reads\n")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.set_dpi(300.0)
    box = ax.get_position()
    ax.set_position([0.1, 0.1, box.width*0.8, box.height])
    min_value, max_value = float("inf"), float("-inf")
    for aln_filename in os.listdir(ALN_EXPL_DIRNAME):
        species_name = aln_filename.replace(".txt", "")
        short_species_name = get_short_name(species_name)
        color = dict_species_color[species_name]

        # Parse alignments to get error rates
        aln_file = open(ALN_EXPL_DIRNAME + aln_filename, "r")
        read_file = open(RAW_READ_DIRNAME + species_name + ".fastq", "r")
        dict_error_quality = {}
        nb_aln_done = 0
        print("  ", species_name)
        while True:
            header = aln_file.readline().rstrip()
            if not header:
                break
            nb_aln_done += 1

            genome = aln_file.readline().rstrip()
            read = aln_file.readline().rstrip()

            read_id = header.split(" ; ")[0].split("=")[1]
            soft_clips = ast.literal_eval(header.split(" ; ")[3].split("=")[1])
            error_rate = get_error_rate(genome, read)
            quality, read_file = get_quality_read(species_name, read_file, read_id, soft_clips)
            if quality not in dict_error_quality:
                dict_error_quality[quality] = [0, 0]
            dict_error_quality[quality][0] += error_rate
            dict_error_quality[quality][1] += 1

            if nb_aln_done > MAX_ALN_NB:
                break

        read_file.close()
        aln_file.close()

        list_qualities, list_error_rates = [], []
        for quality in sorted(dict_error_quality):
            sum_error_rates, occurrences = dict_error_quality[quality]
            if occurrences < MIN_OCC_READ:
                continue
            error_rate = sum_error_rates / occurrences
            list_qualities += [quality]
            list_error_rates += [error_rate]

        plt.plot(list_qualities, list_error_rates, label=short_species_name, color=color)
        if len(list_qualities) == 0:
            print("No data kept. Please consider lower err_qual_min_occ_read value in seqfailr main file.")
            continue
        if min(list_qualities) < min_value:
            min_value = min(list_qualities)
        if max(list_qualities) > max_value:
            max_value = max(list_qualities)

        output_raw.write(f"{species_name}\t{list_qualities}\t{list_error_rates}\n")
    output_raw.close()

    # Reorder legend's labels
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ordered_label_list = [elt for elt in L_ordered_species if elt in by_label]
    ordered_label_values = [by_label[k] for k in ordered_label_list]
    ax.legend(ordered_label_values, ordered_label_list, title="Species", bbox_to_anchor=(1.35,1))

    # Theoretical Phred score: Q = -10 log10(P) ; P = 10**(-Q/10)
    L_val_Q = np.arange(min_value, max_value, 0.5)
    L_val_P = [10**(-q/10) for q in L_val_Q]
    plt.plot(L_val_Q, np.asarray(L_val_P) * 100, "--", color="k")
    plt.xlabel("Mean quality score of reads")
    plt.ylabel("Error rates of reads (%)")
    plt.savefig(OUTPUT_PLOT + "error_rate_quality_score_reads.png")
    plt.close()

def plot_quality_error_for_reads_group():
    """
    Computes error rate depending on quality scores, at read level
    And plot results
    In case species were grouped
    """

    # Get group-species matches
    dict_species_group = {}
    dict_group_color = {}
    group_file = open(FILENAME_SPECIES_GROUPS, "r")
    for line in group_file:
        group_name, species_name, color = line.rstrip().split("\t")
        species_name = species_name.replace(" ", "_")
        dict_species_group[species_name] = group_name
        dict_group_color[group_name] = color
    group_file.close()


    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.set_dpi(300.0)
    box = ax.get_position()
    ax.set_position([0.1, 0.1, box.width*0.8, box.height])
    min_value, max_value = float("inf"), float("-inf")
    dict_to_plot = {}
    for aln_filename in os.listdir(ALN_EXPL_DIRNAME):
        species_name = aln_filename.replace(".txt", "")
        group_name = dict_species_group[species_name.replace(" ", "_")]
        if group_name not in dict_to_plot:
            dict_to_plot[group_name] = {}

        # Parse alignments to get error rates
        aln_file = open(ALN_EXPL_DIRNAME + aln_filename, "r")
        read_file = open(RAW_READ_DIRNAME + species_name + ".fastq", "r")
        dict_error_quality = {}
        nb_aln_done = 0
        print("  ", species_name, group_name)
        while True:
            header = aln_file.readline().rstrip()
            if not header:
                break
            nb_aln_done += 1

            genome = aln_file.readline().rstrip()
            read = aln_file.readline().rstrip()

            read_id = header.split(" ; ")[0].split("=")[1]
            soft_clips = ast.literal_eval(header.split(" ; ")[3].split("=")[1])
            error_rate = get_error_rate(genome, read)
            quality, read_file = get_quality_read(species_name, read_file, read_id, soft_clips)
            if quality not in dict_error_quality:
                dict_error_quality[quality] = [0, 0]
            dict_error_quality[quality][0] += error_rate
            dict_error_quality[quality][1] += 1

            if nb_aln_done > MAX_ALN_NB:
                break

        read_file.close()
        aln_file.close()

        list_qualities, list_error_rates = [], []
        for quality in sorted(dict_error_quality):
            sum_error_rates, occurrences = dict_error_quality[quality]
            if occurrences < MIN_OCC_READ:
                continue
            error_rate = sum_error_rates / occurrences
            if quality not in dict_to_plot[group_name]:
                dict_to_plot[group_name][quality] = []
            dict_to_plot[group_name][quality] += [error_rate]
            list_qualities += [quality]
            list_error_rates += [error_rate]

        if len(list_qualities) == 0:
            print("No data kept. Please consider lower err_qual_min_occ_read value in seqfailr main file.")
            continue
        if min(list_qualities) < min_value:
            min_value = min(list_qualities)
        if max(list_qualities) > max_value:
            max_value = max(list_qualities)

    output_raw = open(OUTPUT_RAW + "error_rate_quality_score_reads.txt", "w")
    output_raw.write("Group_name\tMean quality score reads\tMean error rate reads\n")
    for group_name in dict_to_plot:
        list_qualities, list_error_rates = [], []
        for quality in sorted(dict_to_plot[group_name]):
            list_qualities += [quality]
            list_error_rates += [np.mean(dict_to_plot[group_name][quality])]
        plt.plot(list_qualities, list_error_rates, label=group_name, color=dict_group_color[group_name])
        output_raw.write(f"{group_name}\t{list_qualities}\t{list_error_rates}\n")
    output_raw.close()

    ax.legend(title="Groups", bbox_to_anchor=(1.35,1))


    # Theoretical Phred score: Q = -10 log10(P) ; P = 10**(-Q/10)
    L_val_Q = np.arange(min_value, max_value, 0.5)
    L_val_P = [10**(-q/10) for q in L_val_Q]
    plt.plot(L_val_Q, np.asarray(L_val_P) * 100, "--", color="k")
    plt.xlabel("Mean quality score of reads")
    plt.ylabel("Error rates of reads (%)")
    plt.savefig(OUTPUT_PLOT + "error_rate_quality_score_reads.png")
    plt.close()


def plot_quality_error_for_read_windows():
    """
    Computes error rate depending on quality scores, for WINDOW_LENGTH-bases windows
    And plot results
    In case species were not grouped
    """

    # file to save raw results
    output_raw = open(OUTPUT_RAW + f"error_rate_quality_score_reads_{WINDOW_LENGTH}bases_windows.txt", "w")
    output_raw.write("Species_name\tMean quality score reads\tMean error rate reads\n")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.set_dpi(300.0)
    box = ax.get_position()
    ax.set_position([0.1, 0.1, box.width*0.8, box.height])
    min_value, max_value = float("inf"), float("-inf")

    for aln_filename in os.listdir(ALN_EXPL_DIRNAME):
        species_name = aln_filename.replace(".txt", "")
        print("  ", species_name)
        short_species_name = get_short_name(species_name)
        color = dict_species_color[species_name]
        aln_file = open(ALN_EXPL_DIRNAME + aln_filename, "r")
        read_file = open(RAW_READ_DIRNAME + species_name + ".fastq", "r")
        dict_error_quality = {}
        nb_aln_done = 0
        while True:
            header = aln_file.readline().rstrip()
            if not header:
                break
            nb_aln_done += 1
            if nb_aln_done == MAX_ALN_NB: # limited to MAX_ALN_NB alignments to limit computing time
                break

            genome = aln_file.readline().rstrip() # genome
            read = aln_file.readline().rstrip() # read
            read_id = header.split(" ; ")[0].split("=")[1]
            soft_clips = ast.literal_eval(header.split(" ; ")[3].split("=")[1])
            read_file = compute_err_rate_quality_window(dict_error_quality, species_name,
                                                        genome, read, read_file,
                                                        read_id, soft_clips)
        read_file.close()
        aln_file.close()

        # Plot
        list_qualities, list_error_rates = [], []
        for quality in sorted(dict_error_quality):
            sum_error_rates, occurrences = dict_error_quality[quality]
            if occurrences < MIN_OCC_WINDOW:
                continue
            error_rate = sum_error_rates / occurrences
            list_qualities += [quality]
            list_error_rates += [error_rate]

        plt.plot(list_qualities, list_error_rates, label=short_species_name, color=color)
        if len(list_qualities) == 0:
            continue
        if min(list_qualities) < min_value:
            min_value = min(list_qualities)
        if max(list_qualities) > max_value:
            max_value = max(list_qualities)

        output_raw.write(f"{species_name}\t{list_qualities}\t{list_error_rates}\n")

    output_raw.close()


    # Reorder legend's labels
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ordered_label_list = [elt for elt in L_ordered_species if elt in by_label]
    ordered_label_values = [by_label[k] for k in ordered_label_list]
    ax.legend(ordered_label_values, ordered_label_list, title="Species", bbox_to_anchor=(1.35,1))

    # Theoretical Phred score: Q = -10 log10(P) ; P = 10**(-Q/10)
    L_val_Q = np.arange(min_value, max_value, 0.5)
    L_val_P = [10**(-q/10) for q in L_val_Q]
    plt.plot(L_val_Q, np.asarray(L_val_P) * 100, "--", color="k")
    plt.xlabel("Quality score")
    plt.ylabel("Error rates (%)")

    plt.savefig(OUTPUT_PLOT + f"error_rate_quality_score_reads_{WINDOW_LENGTH}bases_windows.png")

    plt.close()



def plot_quality_error_for_read_windows_group():
    """
    Computes error rate depending on quality scores, for WINDOW_LENGTH-bases windows
    And plot results
    In case species were grouped
    """


    # Get group-species matches
    dict_species_group = {}
    dict_group_color = {}
    group_file = open(FILENAME_SPECIES_GROUPS, "r")
    for line in group_file:
        group_name, species_name, color = line.rstrip().split("\t")
        species_name = species_name.replace(" ", "_")
        dict_species_group[species_name] = group_name
        dict_group_color[group_name] = color
    group_file.close()


    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.set_dpi(300.0)
    box = ax.get_position()
    ax.set_position([0.1, 0.1, box.width*0.8, box.height])
    min_value, max_value = float("inf"), float("-inf")
    dict_to_plot = {}

    for aln_filename in os.listdir(ALN_EXPL_DIRNAME):
        species_name = aln_filename.replace(".txt", "")
        group_name = dict_species_group[species_name.replace(" ", "_")]
        if group_name not in dict_to_plot:
            dict_to_plot[group_name] = {}
        color = dict_group_color[group_name]
        print("  ", species_name, group_name)

        aln_file = open(ALN_EXPL_DIRNAME + aln_filename, "r")
        read_file = open(RAW_READ_DIRNAME + species_name + ".fastq", "r")
        dict_error_quality = {}
        nb_aln_done = 0
        while True:
            header = aln_file.readline().rstrip()
            if not header:
                break
            nb_aln_done += 1
            if nb_aln_done == MAX_ALN_NB: # limited to MAX_ALN_NB alignments to limit computing time
                break

            genome = aln_file.readline().rstrip() # genome
            read = aln_file.readline().rstrip() # read
            read_id = header.split(" ; ")[0].split("=")[1]
            soft_clips = ast.literal_eval(header.split(" ; ")[3].split("=")[1])
            read_file = compute_err_rate_quality_window(dict_error_quality, species_name,
                                                        genome, read, read_file,
                                                        read_id, soft_clips)
        read_file.close()
        aln_file.close()

        # Plot
        list_qualities, list_error_rates = [], []
        for quality in sorted(dict_error_quality):
            sum_error_rates, occurrences = dict_error_quality[quality]
            if occurrences < MIN_OCC_WINDOW:
                continue
            error_rate = sum_error_rates / occurrences
            if quality not in dict_to_plot[group_name]:
                dict_to_plot[group_name][quality] = []
            dict_to_plot[group_name][quality] += [error_rate]
            list_qualities += [quality]
            list_error_rates += [error_rate]

        if len(list_qualities) == 0:
            continue
        if min(list_qualities) < min_value:
            min_value = min(list_qualities)
        if max(list_qualities) > max_value:
            max_value = max(list_qualities)


    output_raw = open(OUTPUT_RAW + f"error_rate_quality_score_reads_{WINDOW_LENGTH}bases_windows.txt", "w")
    output_raw.write("Group_name\tMean quality score reads\tMean error rate reads\n")
    for group_name in dict_to_plot:
        list_qualities, list_error_rates = [], []
        for quality in sorted(dict_to_plot[group_name]):
            list_qualities += [quality]
            list_error_rates += [np.mean(dict_to_plot[group_name][quality])]
        plt.plot(list_qualities, list_error_rates, label=group_name, color=dict_group_color[group_name])
        output_raw.write(f"{group_name}\t{list_qualities}\t{list_error_rates}\n")
    output_raw.close()


    ax.legend(title="Species", bbox_to_anchor=(1.35,1))

    # Theoretical Phred score: Q = -10 log10(P) ; P = 10**(-Q/10)
    L_val_Q = np.arange(min_value, max_value, 0.5)
    L_val_P = [10**(-q/10) for q in L_val_Q]
    plt.plot(L_val_Q, np.asarray(L_val_P) * 100, "--", color="k")
    plt.xlabel("Quality score")
    plt.ylabel("Error rates (%)")

    plt.savefig(OUTPUT_PLOT + f"error_rate_quality_score_reads_{WINDOW_LENGTH}bases_windows.png")

    plt.close()


# --------------------------------------------------------------------------------------------------
# Main

if __name__ == "__main__":

    # Parses arguments
    NUMBER_EXPECTED_ARGUMENTS = 10
    if len(sys.argv) != NUMBER_EXPECTED_ARGUMENTS + 1:
        print(f"   ERROR: Wrong number of arguments: {NUMBER_EXPECTED_ARGUMENTS} expected but {len(sys.argv)-1} given.")
        sys.exit(2)
    ALN_EXPL_DIRNAME, RAW_READ_DIRNAME, OUTPUT_RAW, OUTPUT_PLOT, FILENAME_SPECIES_GC_COLOR, FILENAME_SPECIES_GROUPS, WINDOW_LENGTH, MAX_ALN_NB, MIN_OCC_READ, MIN_OCC_WINDOW = sys.argv[1:]
    WINDOW_LENGTH = int(WINDOW_LENGTH)
    MAX_ALN_NB = int(MAX_ALN_NB)
    MIN_OCC_READ = int(MIN_OCC_READ)
    MIN_OCC_WINDOW = int(MIN_OCC_WINDOW)

    if MAX_ALN_NB == -1:
        MAX_ALN_NB = float("inf")


    # Get color associated to each species, and list of ordered species (according to GC content)
    dict_species_color, L_ordered_species = get_color(FILENAME_SPECIES_GC_COLOR)


    # --- 1/ Error rate depending on quality score, at read level ---
    print("   - Computes error rate depending on quality scores, at read level")
    print(f"      A minimum occurrence number ({MIN_OCC_READ}) is required")
    if os.path.exists(FILENAME_SPECIES_GROUPS):
        plot_quality_error_for_reads_group()
    else:
        plot_quality_error_for_reads()


    # --- 2/ Error rate depending on quality score, for WINDOW_LENGTH-bases windows scale
    # WARNING: limited to max MAX_ALN_NB alignments per species (considered enough) to limit computing time
    print()
    print(f"   - Computes error rate depending on quality scores, at {WINDOW_LENGTH}-bases window scale")
    print(f"      A minimum occurrence number ({MIN_OCC_WINDOW}) is required")
    print(f"     Warning: limit to max {MAX_ALN_NB} alignments per species = enough to get relevant data" +
          " and limiting running time of the script")
    if os.path.exists(FILENAME_SPECIES_GROUPS):
        plot_quality_error_for_read_windows_group()
    else:
        plot_quality_error_for_read_windows()

