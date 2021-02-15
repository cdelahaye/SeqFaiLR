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

params = {'legend.fontsize': 16,
          'legend.title_fontsize': 16,
          'legend.labelspacing': 0.1,
          'legend.borderpad': 0.3,
          'figure.figsize': (14, 9),
          'axes.labelsize': 16,
          'axes.titlesize': 16,
          'xtick.labelsize': 16,
          'ytick.labelsize': 16}
plt.rcParams.update(params)


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

def get_error_rate():
    """ Computes error rate for the current alignment
    """
    error_count = 0
    for i, base_genome in enumerate(genome):
        base_read = read[i]
        if base_read != base_genome:
            error_count += 1
    error_rate = round(error_count / len(genome) * 100, 4)
    return error_rate


def get_quality_read(read_file):
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


def compute_err_rate_quality_window(read_file):
    """ Computes error rates and gc content for current alignment,
    for 100-bases sliding windows along read
    Updates dictionary storing results
    """
    WINDOW_LENGTH = 100

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
        if quality_score_mean not in dict_error_quality:
            dict_error_quality[quality_score_mean] = [0, 0]
        dict_error_quality[quality_score_mean][0] += error_rate
        dict_error_quality[quality_score_mean][1] += 1

        # Update quality portion and score
        quality_score = quality_score - (ord(quality_str[0]) - 33)
        quality_str = quality_str[1:]
        quality_str_part = quality_str[0: WINDOW_LENGTH]
        quality_score = quality_score + (ord(quality_str_part[-1]) - 33)

        # Update start positions
        start_pos += 1
        while read[start_pos] == "-":
            start_pos += 1


# --------------------------------------------------------------------------------------------------
# Main

if __name__ == "__main__":

    ## Parses arguments
    if len(sys.argv) != 6:
        print(f"ERROR: Wrong number of arguments: 5 expected but {len(sys.argv)-1} given.")
        sys.exit(2)
    ALN_EXPL_DIRNAME = sys.argv[1]
    RAW_READ_DIRNAME = sys.argv[2]
    OUTPUT_RAW = sys.argv[3]
    OUTPUT_PLOT = sys.argv[4]
    FILENAME_SPECIES_GC = sys.argv[5]

    if ALN_EXPL_DIRNAME[-1] != "/":
        ALN_EXPL_DIRNAME += "/"
    if RAW_READ_DIRNAME[-1] != "/":
        RAW_READ_DIRNAME += "/"
    if OUTPUT_PLOT[-1] != "/":
        OUTPUT_PLOT += "/"
    if OUTPUT_RAW[-1] != "/":
        OUTPUT_RAW += "/"


    # Get color for each species and a list of ordered species (depending on GC content)
    dict_species_color = {}
    dict_gc_species = {}
    file = open(FILENAME_SPECIES_GC, "r")
    for line in file:
        species_name, color, gc = line.rstrip().split(" ; ")
        gc = float(gc)
        species_name = species_name.replace(" ", "_")
        if color[0] != "#":
            color = "#" + color
        dict_species_color[species_name] = color
        if gc not in dict_gc_species:
            dict_gc_species[gc] = []
        dict_gc_species[gc] += [species_name]
    file.close()
    list_ordered_species_gc = []
    for gc in sorted(dict_gc_species):
        for sp in dict_gc_species[gc]:
            list_ordered_species_gc += [get_short_name(sp)]


    # First graph: error rate depending on quality score, at read level
    print("Computes error rate depending on quality scores, at read level")

    # file to save raw results
    output_raw = open(OUTPUT_RAW + "error_rate_quality_score_reads.txt", "w")
    output_raw.write("Species_name\tMean quality score reads\tMean error rate reads\n")

    fig, ax = plt.subplots()
    fig.set_dpi(300.0)
    for aln_filename in os.listdir(ALN_EXPL_DIRNAME):
        species_name = aln_filename.replace(".txt", "")
        short_species_name = get_short_name(species_name)
        color = dict_species_color[species_name]
        aln_file = open(ALN_EXPL_DIRNAME + aln_filename, "r")
        read_file = open(RAW_READ_DIRNAME + species_name + ".fastq", "r")
        dict_error_quality = {}
        nb_aln_done = 0
        print(species_name)
        while True:
            header = aln_file.readline().rstrip()
            if not header:
                break
            genome = aln_file.readline().rstrip()
            read = aln_file.readline().rstrip()

            nb_aln_done += 1

            # Bound number of alignment analysed
            if nb_aln_done > 100000:
                break

            read_id = header.split(" ; ")[0].split("=")[1]
            if "reverse" in header:
                strand = "reverse"
            else:
                strand = "forward"
            soft_clips = ast.literal_eval(header.split(" ; ")[-1].split("=")[1])
            error_rate = get_error_rate()
            quality, read_file = get_quality_read(read_file)
            if quality not in dict_error_quality:
                dict_error_quality[quality] = [0, 0]
            dict_error_quality[quality][0] += error_rate
            dict_error_quality[quality][1] += 1
        read_file.close()
        aln_file.close()

        list_qualities, list_error_rates = [], []
        if "Human" in aln_filename:
            threshold = 10000
        else:
            threshold = 10
        for quality in sorted(dict_error_quality):
            sum_error_rates, occurrences = dict_error_quality[quality]
            if occurrences < threshold:
                continue
            error_rate = sum_error_rates / occurrences
            list_qualities += [quality]
            list_error_rates += [error_rate]

        plt.plot(list_qualities, list_error_rates, label=short_species_name, color=color)
        output_raw.write(f"{species_name}\t{list_qualities}\t{list_error_rates}\n")
    output_raw.close()

    # Reorder legend's labels
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ordered_label_list = [elt for elt in list_ordered_species_gc if elt in by_label]
    ordered_label_values = [by_label[k] for k in ordered_label_list]
    plt.legend(ordered_label_values, ordered_label_list, title="Species")

    # Theoretical Phred score: Q = -10 log10(P)
    L_val_P = np.arange(0.05, 12, 0.5)/100
    L_val_Q = -10 * np.log10(L_val_P)
    plt.plot(L_val_Q, np.asarray(L_val_P) * 100, "--", color="k")
    plt.text(14, 1, 'Expected Phred \nscores distribution', style='italic',
             ha="center", fontsize=18)
    plt.xlim(10, 35)
    plt.xticks(ticks=np.arange(10, 35, 2), labels=np.arange(10, 41, 2))
    plt.ylim(0, 12)
    plt.yticks(ticks=np.arange(0, 12, 1), labels=np.arange(0, 12, 1))
    plt.xlabel("Mean quality score of reads")
    plt.ylabel("Error rates of reads (%)")
    plt.savefig(OUTPUT_PLOT + "error_rate_quality_score_reads.png")
    plt.close()








    # Second graph: error rate depending on quality score, for 100-base windows scale
    # WARNING: limited to max 10,000 alignments per species (considered enough)
    print("Computes error rate depending on quality scores, at 100-bases window scale")
    print("Warning: limit to max 10,000 alignments per species = enough to get relevant data" +
          " and limiting running time of the script")
    # file to save raw results
    output_raw = open(OUTPUT_RAW + "error_rate_quality_score_reads_100bases_windows.txt", "w")
    output_raw.write("Species_name\tMean quality score reads\tMean error rate reads\n")

    fig, ax = plt.subplots()
    fig.set_dpi(300.0)
    for aln_filename in os.listdir(ALN_EXPL_DIRNAME):
        species_name = aln_filename.replace(".txt", "")
        short_species_name = get_short_name(species_name)
        color = dict_species_color[species_name]
        aln_file = open(ALN_EXPL_DIRNAME + aln_filename, "r")
        read_file = open(RAW_READ_DIRNAME + species_name + ".fastq", "r")
        dict_error_quality = {}
        nb_aln_done = 0
        print(species_name)
        while True:
            header = aln_file.readline().rstrip()
            if not header:
                break
            genome = aln_file.readline().rstrip()
            read = aln_file.readline().rstrip()

            nb_aln_done += 1

            if nb_aln_done % 500 == 0:
                print(f"{nb_aln_done} alignments processed")
            if nb_aln_done == 10000:
                break

            read_id = header.split(" ; ")[0].split("=")[1]
            if "reverse" in header:
                strand = "reverse"
            else:
                strand = "forward"
            soft_clips = ast.literal_eval(header.split(" ; ")[-1].split("=")[1])

            read_file = compute_err_rate_quality_window(read_file)

        read_file.close()
        aln_file.close()

        # Plot
        list_qualities, list_error_rates = [], []
        if "Human" in aln_filename:
            threshold = 100000
        else:
            threshold = 1000
        for quality in sorted(dict_error_quality):
            sum_error_rates, occurrences = dict_error_quality[quality]
            if occurrences < threshold:
                continue
            error_rate = sum_error_rates / occurrences
            list_qualities += [quality]
            list_error_rates += [error_rate]

        plt.plot(list_qualities, list_error_rates, label=short_species_name, color=color)
        # Reorder legend's labels
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ordered_label_list = [elt for elt in list_ordered_species_gc if elt in by_label]
        ordered_label_values = [by_label[k] for k in ordered_label_list]
        plt.legend(ordered_label_values, ordered_label_list, title="Species")

        # Theoretical Phred score: Q = -10 log10(P)
        L_val_P = np.arange(0.01, 42, 0.5)/100
        L_val_Q = -10 * np.log10(L_val_P)
        plt.plot(L_val_Q, np.asarray(L_val_P) * 100, "--", color="k")
        plt.text(9, 1, 'Expected Phred \nscores distribution', style='italic',
                 ha="center", fontsize=18)
        plt.xlim(0, 40)
        plt.xticks(ticks=np.arange(0, 41, 5), labels=np.arange(0, 41, 5))
        #plt.ylim(0, 35)
        #plt.yticks(ticks=np.arange(0, 36, 5), labels=np.arange(0, 36, 5))
        plt.xlabel("Quality score")
        plt.ylabel("Error rates (%)")

        plt.savefig(OUTPUT_PLOT + "error_rate_quality_score_reads_100bases_windows.png")
        output_raw.write(f"{species_name}\t{list_qualities}\t{list_error_rates}\n")
    output_raw.close()
    plt.close()
