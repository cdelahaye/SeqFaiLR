#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Analyses errors in homopolymeric regions, and computes (and plots) several information:
    - quantifying homopolymers in reference genomes
    - quantifying errors related to homopolymeric regions (compared to all errors)
    - error rates (mismatch, insertion, deletion) depending on homopolymer length
    - differences between expected and sequenced homopolymer lengths
All plots are made here by distinguishing groups (as defined by user)
"""

# --------------------------------------------------------------------------------------------------
# Packages

import sys
import os
import re
import time
import datetime
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
# Parameters for output plots
PARAMS = {'legend.fontsize': 16,
          'legend.title_fontsize': 16,
          'legend.labelspacing': 0.1,
          'legend.borderpad': 0.3,
          'legend.columnspacing': 1,
          'legend.handletextpad': 0.5,
          'legend.handlelength': 0.8,
          'figure.figsize': (10, 8),
          'axes.labelsize': 16,
          'axes.titlesize': 16,
          'xtick.labelsize': 16,
          'ytick.labelsize': 16}
plt.rcParams.update(PARAMS)

# --------------------------------------------------------------------------------------------------
# Functions

def get_groups(filename):
    """
    Returns a dictionary of species names (keys), associated with their group name (value)
    Also a similar dictionary that associated species names (keys) with color of their group (value)
    """
    dictionary_group = {}
    dictionary_color = {}
    file = open(filename, "r")
    for line in file:
        group_name, species_name, color = line.rstrip().split("\t")
        species_name = species_name.split(".")[0]
        dictionary_group[species_name] = group_name
        dictionary_color[group_name] = color
    file.close()

    return dictionary_group, dictionary_color

def get_short_name(long_name):
    L_name = long_name.replace("_", " ").split(" ")
    L_name[0] = L_name[0][0] + "."
    short_name = " ".join(L_name)
    return short_name



def build_group_per_species(filename):
    """
    If no group were defined by user, create N groups (1 for each of the N species)
    Also returns list of ordered species (based on their GC content)
    """
    dictionary_group = {}
    dictionary_color = {}
    dict_gc = {}
    file = open(filename, "r")
    for line in file:
        species_name, color, gc = line.rstrip().split("\t")
        gc = float(gc)
        species_name = get_short_name(species_name)
        group_name = species_name
        dictionary_group[species_name] = group_name
        dictionary_color[group_name] = color
        if gc not in dict_gc:
            dict_gc[gc] = []
        dict_gc[gc] += [species_name]
        
    L_ordered_species = [get_short_name(species_name) for gc in sorted(dict_gc)
                                                      for species_name in dict_gc[gc]]
    file.close()

    return dictionary_group, dictionary_color, L_ordered_species


def get_total_number_of_lines(filename: str) -> int:
    """Returns number of lines of filename given as input
    """
    with open(filename, "r") as file:
        return sum(1 for _ in file)


def display_progressing_bar(prct: int, current_time: float):
    """ Display a simple progressing of the script, in STDOUT

    Parameters
    ----------
    prct: int
        Percentage of work done

    current_time: float
        Current time at function's calling
        Used to compute elapsed time, and estimate remaining time
    """

    # Percentage of work done and remaining
    prct_done = int(prct//2)
    prct_remaining = int((100/2) - prct_done)

    # Progressing bar
    progressing_bar = "|" + (prct_done)*"#" + (prct_remaining)*"-" + "|"

    # Estimation of remaining time to complete
    elapsed_time = current_time - STARTING_TIME
    estimated_remaining_time = (NB_TOT_ALN-nb_aln_done) * elapsed_time / nb_aln_done
    estimated_remaining_time = str(datetime.timedelta(seconds=round(estimated_remaining_time)))

    # Write to stdout
    sys.stdout.write(f"\r{progressing_bar} ({prct}% ; {estimated_remaining_time})")
    sys.stdout.flush()


def initiate_dictionary_error_length() -> dict:
    """ Initiate a dictionary that will store errors depending on their type (mismatch, insertion,
        deletion) and their length, for each species group
        Mismatches and indels are not stored explicitely, but in the following order:
            mismatch (position 0), insertion (position 1) and deletion (position 2)
    """
    dictionary = {}
    for group_name in L_groups:
        dictionary[group_name] = {}
        for error_length in L_lengths:
            dictionary[group_name][error_length] = [0, 0, 0] # mismatch, insertion, deletion
    return dictionary

def initiate_diff_len_dict() -> dict:
    """Add a 'sub-dictionary' for the current species, which will, for each species group and
    for each genomic homopolymer length, store the sequenced homopolymer length (i.e. in the read)
    """
    dictionary = {}
    for group_name in L_groups:
        dictionary[group_name] = {}
        for error_length in L_lengths:
            dictionary[group_name][error_length] = {}
    return dictionary


def initiate_prct_correct_dict():
    """

    """
    dictionary = {}
    dictionary_detail = {}
    for group_name in L_groups:
        dictionary[group_name] = {}
        dictionary_detail[group_name] = {}
        for length in L_lengths:
            dictionary[group_name][length] = []
            for base in L_bases:
                if base not in dictionary_detail[group_name]:
                    dictionary_detail[group_name][base] = {}
                dictionary_detail[group_name][base][length] = []

    return dictionary, dictionary_detail




def get_genomic_homopolymer_distribution(filename, group_name):
    """Compute homopolymer distribution in the current reference genome, and
    updates the associated dictionary
    """

    # --- Add "space" to store results ---
    for length in dict_genomic_distrib[group_name]:
        for base in L_bases:
            if base not in dict_genomic_distrib[group_name][length]:
                dict_genomic_distrib[group_name][length][base] = [0]
            else:
                dict_genomic_distrib[group_name][length][base] += [0]


    # --- Compute homopolymer distribution ---
    with open(filename, "r") as reference_genome_file:
        genome = ""
        reference_genome_file.readline() # skip first line = header
        while True:
            line = reference_genome_file.readline().rstrip()
            if not line:
                break
            if line[0] == ">": # Compute results when new chromosome
                for res in re.finditer(pattern_homopolymer, genome):
                    start, end = res.span()
                    homopolymer_length = end - start
                    if homopolymer_length >= MAX_HOMOPOLYMER_LENGTH:
                        homopolymer_length = f"{MAX_HOMOPOLYMER_LENGTH}+"
                    base = res.group()[0]
                    dict_genomic_distrib[group_name][homopolymer_length][base][-1] += 1
                genome = ""
                continue
            genome += line.upper()

        # Compute last results
        for res in re.finditer(pattern_homopolymer, genome):
            start, end = res.span()
            homopolymer_length = end - start
            if homopolymer_length >= MAX_HOMOPOLYMER_LENGTH:
                homopolymer_length = f"{MAX_HOMOPOLYMER_LENGTH}+"
            base = res.group()[0]
            dict_genomic_distrib[group_name][homopolymer_length][base][-1] += 1


    # --- Plot ---
    fig, ax = plt.subplots()
    fig.set_dpi(300.0)
    bar_width = 0.42

    dict_to_plot = {}
    for group_name in dict_genomic_distrib:
        if group_name not in dict_to_plot:
            dict_to_plot[group_name] = {}
        for length in dict_genomic_distrib[group_name]:
            if length not in dict_to_plot[group_name]:
                dict_to_plot[group_name][length] = {}
            for base in dict_genomic_distrib[group_name][length]:
                if base not in dict_to_plot[group_name][length]:
                    dict_to_plot[group_name][length][base] = []
                dict_to_plot[group_name][length][base] += dict_genomic_distrib[group_name][length][base]

    for group_name in dict_to_plot:
        if dict_to_plot[group_name][MIN_HOMOPOLYMER_LENGTH] == {}:
            continue
        color = dict_species_group_color[group_name]
        L_occurences = []
        for length in L_lengths:
            L_occ_group = []
            for base in L_bases:
                if base not in dict_to_plot[group_name][length]:
                    L_occ_group += [0]
                    continue
                median_occurrences = np.median(dict_to_plot[group_name][length][base])
                L_occ_group += [median_occurrences]
            L_occurences += [L_occ_group]

        offset = -0.75 # margin between two bases
        step = 3 # margin between two lengths
        for i, base in enumerate(L_bases):
            color = dict_bases_colors[base]
            L_to_plot_y = [elt[i] for elt in L_occurences]
            plt.bar(np.arange(0, len(L_lengths))*step+offset, L_to_plot_y,
                    width=bar_width, color=color, label=base)
            offset += 0.5

        ax.set_yscale('log')
        ax.set_xlabel("Homopolymer length")
        ax.set_ylabel("Number of occurrences")
        L_xticks = list(np.arange(0, len(L_lengths)*step, step))
        L_xtick_labels = list(np.arange(MIN_HOMOPOLYMER_LENGTH, MAX_HOMOPOLYMER_LENGTH)) \
                            + [f"{MAX_HOMOPOLYMER_LENGTH}+"]
        plt.xticks(L_xticks, L_xtick_labels)
        plt.legend(title="Bases", ncol=len(L_bases))
        plt.title(f"Distribution of homopolymer lengths for {group_name}")
        plt.savefig(OUTPUT_PLOT + f"homopolymer_genomic_distribution_{group_name.replace(' ', '_')}.png")
        plt.close()


    # --- Write raw results in .txt file ---
    RAW_OUTPUT_FILE = open(OUTPUT_RAW + "homopolymer_genomic_distribution.txt", "w")
    for group_name in dict_genomic_distrib:
        RAW_OUTPUT_FILE.write(f"{group_name}\n")
        RAW_OUTPUT_FILE.write("\t".join(["Homopolymer length"] + L_bases) + "\n")
        for homopolymer_length in dict_genomic_distrib[group_name]:
            if dict_genomic_distrib[group_name][homopolymer_length] == {}:
                L_values = ["-" for base in L_bases]
            else:
                L_values = [np.median(dict_genomic_distrib[group_name][homopolymer_length][base]) for base in L_bases]
            RAW_OUTPUT_FILE.write(str(homopolymer_length) + "\t" + "\t".join(map(str, L_values)) + "\n")
    RAW_OUTPUT_FILE.close()



def update_errors(genome_aln, read_aln, group_name, category):
    """
    Computes errors in alignment of genome_aln versus read_aln
    Updates dict_quantify_errors for the associated group_name
    Category can be "All_errors" or "Homopolymer_errors"

    For the N-bases in the genome, we consider the aligned bases as non-erroneous as the
      reference is unknown
    """
    for i, genome_base in enumerate(genome_aln):
        read_base = read_aln[i]
        if genome_base == "N" or genome_base == read_base:
            continue
        if genome_base == "-":
            error_type = "Insertion"
        elif read_base == "-":
            error_type = "Deletion"
        elif genome_base != read_base:
            error_type = "Mismatch"
        dict_quantify_errors[group_name][category][error_type] += 1
        dict_quantify_errors[group_name][category]["Total"] += 1



def update_errors_homopolymer_length(g: str, r: str, homopolymer_length: int):
    """Computes errors in alignment of g (genome) versus r (read)
    Updates global dictionary that store these results, for length homopolymer_length
    Mismatches and indels are not stored explicitely, but in the following order:
        mismatch (position 0), insertion (position 1) and deletion (position 2)
    """
    for i, genome_base in enumerate(g):
        read_base = r[i]
        if genome_base == "N" or genome_base == read_base:
            continue
        if genome_base == "-":
            dict_error_length[group_name][homopolymer_length][1] += 1
        elif read_base == "-":
            dict_error_length[group_name][homopolymer_length][2] += 1
        elif genome_base != read_base:
            dict_error_length[group_name][homopolymer_length][0] += 1


def update_prct_correct_homopolymer(g: str, r: str, group_name:str, homopolymer_length: int):
    """Updates dictionary prct_correct_homopolymer_dict with counts of correctly sequenced
    homopolymers, and total number of homopolymers
    """
    if g == r:
        dict_prct_correct[group_name][homopolymer_length][-1][0] += 1
    dict_prct_correct[group_name][homopolymer_length][-1][1] += 1

def update_prct_correct_homopolymer_detail(g: str, r: str, group_name: str, homopolymer_length: int,
                                           homopolymer_base: str):
    """Same as previous, but for the detailed version
    """
    if g == r:
        dict_prct_correct_detail_bases[group_name][homopolymer_base][homopolymer_length][-1][0] += 1
    dict_prct_correct_detail_bases[group_name][homopolymer_base][homopolymer_length][-1][1] += 1


def update_diff_len_dict(g: str, r: str, group_name: str, homopolymer_length_genome: int,
                         homopolymer_base: str):
    """Updates dictionary that stores homopolymer length differences (genomic and sequenced)
    """
    homopolymer_length_read = len(max(re.compile(homopolymer_base + "*").findall(r)))
    if homopolymer_length_read not in dict_diff_length_sequenced[group_name][homopolymer_length_genome]:
        dict_diff_length_sequenced[group_name][homopolymer_length_genome][homopolymer_length_read] = 0
    dict_diff_length_sequenced[group_name][homopolymer_length_genome][homopolymer_length_read] += 1


def get_statistics(dictionary):
    """Computes minimum, maximum, median and quartiles 1 and 3
    from a given dictionary of the form:
        {expected_value1: {sequenced_value1: occurrences,
                          sequenced_value2: occurrences,
                          ...},
        expected_value2: {sequenced_value1: occurrences,
                          sequenced_value2: occurrences,
                          ...},
        ...}
    Better memory usage than storing a list of occurrences for each expected value
    """
    if dictionary == {}:
        return({'med': -1, 'q1': -1, 'q3': -1, 'whislo': -1, 'whishi': -1})
    # Get N = total number of "observation"
    N = 0
    for homopol_read in dictionary:
        N += dictionary[homopol_read]

    # Minimum, maximum, median and quartiles:
    minimum = sorted(list(dictionary))[0]
    maximum = sorted(list(dictionary))[-1]
    N_q1 = N / 4
    N_med = N / 2
    N_q3 = N * 3/4
    q1, median, q3 = None, None, None
    n0 = 0
    n1 = n0
    for homopol_read in sorted(dictionary):
        n0 = n1
        n1 = n0 + dictionary[homopol_read]
        # First quartile
        if q1 == None:
            if n0 < N_q1 < n1:
                q1 = homopol_read
            elif N_q1 < n1:
                if n0 == N_q1:
                    q1 = homopol_read
        # Median
        if median is None:
            if n0 < N_med < n1:
                median = homopol_read
            elif N_med < n1:
                if n0 == N_med:
                    median = homopol_read
        # Third quartile
        if q3 is None:
            if n0 < N_q3 < n1:
                q3 = homopol_read
            elif N_q3 < n1:
                if n0 == N_q3:
                    q3 = homopol_read
    # Turns 0 values to 1 for log-scaled output graph
    if minimum == 0:
        minimum = 0.9
    return({'med': median, 'q1': q1, 'q3': q3, 'whislo': minimum, 'whishi': maximum})


def try_extend_genomic_homopolymer(start_pos, end_pos, sequence):
    """ If sequence[start: end] start/end by '-' or homopolymer base, extend the
    start and end positions (delimiting genomic homopolymer)
    """
    if sequence[start_pos] in ("-", base_homopolymer):
        while start_pos-1 >= 0 and sequence[start_pos-1] in ("-", base_homopolymer):
            start_pos -= 1
    if sequence[end_pos-1] in ("-", base_homopolymer):
        while end_pos < len(sequence) and sequence[end_pos] in ("-", base_homopolymer):
            end_pos += 1
    while sequence[start_pos] == "-":
        start_pos += 1
    while sequence[end_pos-1] == "-":
        end_pos -= 1
    return start_pos, end_pos


def try_extend_read_homopolymer(start_pos, end_pos, sequence):
    """ If sequence[start: end] start/end by '-' or homopolymer base, extend the
    start and end positions (delimiting read homopolymer)
    """
    if sequence[start_pos: end_pos].replace("-", "") == "":
        return start_pos, end_pos

    if sequence[start_pos] == base_homopolymer:
        while start_pos-1 >= 0 and sequence[start_pos-1] == base_homopolymer:
            start_pos -= 1
    if sequence[end_pos-1] == base_homopolymer:
        while end_pos < len(sequence) and sequence[end_pos] == base_homopolymer:
            end_pos += 1

    return start_pos, end_pos


def compute_results():
    """Computes all results, and output as graph and raw file
    """

    # 1/ --- Global quantification of homopolymer errors. Outputs raw .txt file ---
    RAW_OUTPUT_FILE = open(OUTPUT_RAW + "quantify_errors_homopolymer_vs_total.txt", "w")

    header = "\t".join(L_errors)
    for group_name in L_groups:
        RAW_OUTPUT_FILE.write(f"\n{group_name}\n")
        L_homopolymer_errors, L_all_errors, L_ratio_homopolymer_errors = [], [], []
        # gather results for each error type
        for error in L_errors:
            nb_homopolymer_errors = round(dict_quantify_errors[group_name]["Homopolymer_errors"][error], 2)
            nb_all_errors = round(dict_quantify_errors[group_name]["All_errors"][error], 2)
            if nb_all_errors == 0:
                ratio_homopolymer_errors = 0
            else:
                ratio_homopolymer_errors = round(nb_homopolymer_errors / nb_all_errors * 100, 2)
            L_homopolymer_errors += [nb_homopolymer_errors]
            L_all_errors += [nb_all_errors]
            L_ratio_homopolymer_errors += [ratio_homopolymer_errors]
        # write in output file
        homopolymer_errors_str = "\t".join(["H: Homopol."] + \
                                            list(map(str, L_homopolymer_errors)) + \
                                            [str(round(np.median(L_homopolymer_errors), 2))])
        all_errors_str = "\t".join(["A: All errors"] + \
                                   list(map(str, L_all_errors)) + \
                                   [str(round(np.median(L_all_errors), 2))])
        ratio_homopolymer_errors = "\t".join(["Ratio H/A (%)"] + \
                                              list(map(str, L_ratio_homopolymer_errors)) + \
                                              [str(round(np.median(L_ratio_homopolymer_errors), 2))])
        RAW_OUTPUT_FILE.write("\n".join([header, homopolymer_errors_str,
                                         all_errors_str, ratio_homopolymer_errors]) + "\n")
    RAW_OUTPUT_FILE.close()



    # 2/ --- Computes a pie of error rates (mismatches and indel) depending on homopolymer length ---
    #   a) Plot for each group
    L_color = [COLOR_MISMATCH, COLOR_INSERTION, COLOR_DELETION]
    L_label = ["Mismatches", "Insertions", "Deletions"]
    nb_max_row = 3
    nb_max_column = (len(L_lengths) // nb_max_row) + 1
    for group_name in L_groups:
        i_row, i_column = 0, 0
        fig = plt.figure()
        for length in dict_error_length[group_name]:

            ax = plt.subplot2grid((nb_max_row, nb_max_column), (i_row, i_column))
            wedges, _, _ = ax.pie(dict_error_length[group_name][length],
                                  colors=L_color, autopct='%1.0f%%',
                                  pctdistance=0.8,
                                  textprops=dict(color="k", size=12),
                                  radius=1.2)
            title = "Length " + str(length)
            plt.xlabel(title)
            if i_column < nb_max_column - 1:
                i_column += 1
            else:
                i_row += 1
                i_column = 0

        fig.legend(title="Error:", labels=L_label, ncol=len(L_label))
        plt.savefig(OUTPUT_PLOT + f"error_rates_depending_on_homopolymer_length_{group_name.replace(' ', '_')}.png")
        plt.close()

    #   b) Save raw results in .txt file
    RAW_OUTPUT_FILE = open(OUTPUT_RAW + "error_rates_depending_on_homopolymer_length.txt", "w")
    for group_name in L_groups:
        RAW_OUTPUT_FILE.write(f"\n{group_name}\n")
        L_error = ["Mismatch", "Insertion", "Deletion"]
        RAW_OUTPUT_FILE.write("\t".join(["Homopolymer length"] + L_error) + "\n")
        for homopolymer_length in dict_error_length[group_name]:
            RAW_OUTPUT_FILE.write(str(homopolymer_length) + "\t" +
                                  "\t".join(map(str, dict_error_length[group_name][homopolymer_length]))
                                  + "\n")
    RAW_OUTPUT_FILE.close()



    # 3/ --- Plots ratio of correctly sequenced homopolymers over and total number of homopolymers ---
    # -- if user defined groups: boxplot --
    if os.path.exists(FILE_SPECIES_GROUP):
        # Prepare data for the plot, and converts counts to ratios
        nb_groups = len(L_groups)
        offset_intra = 0.25 # margin between groups (same length)
        offset_inter = (nb_groups+2) * offset_intra  # margin between homopolymer length
        position = 0
        L_patches = []
        for i_grp, group_name in enumerate(dict_prct_correct):
            position = i_grp * offset_intra
            color = dict_species_group_color[group_name]
            L_positions = [] # Position of plot
            L_prct_correct = [] # Values to plot
            for length in dict_prct_correct[group_name]:
                L_positions += [position]
                position += offset_inter
                tmp_list = []
                for i, _ in enumerate(dict_prct_correct[group_name][length]):
                    if dict_prct_correct[group_name][length][i][1] != 0:
                        tmp_list += [(dict_prct_correct[group_name][length][i][0] /
                                      dict_prct_correct[group_name][length][i][1] * 100)]
                L_prct_correct += [tmp_list]

            box = plt.boxplot(L_prct_correct, positions=L_positions,
                              widths=BOXPLOT_WIDTH, patch_artist=True)
            for item in ['medians']:
                plt.setp(box[item], color="black")
            for patch in box['boxes']:
                patch.set(facecolor=color)
            patch = mpatches.Patch(color=color, label=group_name)
            L_patches += [patch]
        plt.xticks(ticks=[i* offset_inter + ((nb_groups-1)/2 * offset_intra) for i in range(len(L_lengths))],
                   labels=L_lengths)
        plt.xlabel("Theoretical homopolymer length (reference genome)")
        plt.ylabel("Errorless sequenced homopolymer ratio (reads)")
        plt.legend(title="Species:", handles=L_patches)

    # -- otherwise user defined groups: plot --
    else:
        # Prepare data for the plot, and converts counts to ratios
        nb_groups = len(L_groups)
        offset_intra = 0.401 # margin between groups (same length)
        offset_inter = (nb_groups+2) * offset_intra  # margin between homopolymer length
        position = 0
        for i_grp, group_name in enumerate(L_ordered_species):
            position = i_grp * offset_intra
            color = dict_species_group_color[group_name]
            L_positions = [] # Position of plot
            L_prct_correct = [] # Values to plot
            for length in dict_prct_correct[group_name]:
                L_positions += [position]
                position += offset_inter
                if dict_prct_correct[group_name][length] == []:
                    tmp_list = [float("-inf")]
                else:
                    tmp_list = []
                    for i, _ in enumerate(dict_prct_correct[group_name][length]):
                        if dict_prct_correct[group_name][length][i][1] != 0:
                            tmp_list += [(dict_prct_correct[group_name][length][i][0] /
                                          dict_prct_correct[group_name][length][i][1] * 100)]
                        else:
                            tmp_list += [float("-inf")]
                L_prct_correct += [tmp_list]
            plt.plot(L_positions, L_prct_correct, "o", color=color, label=group_name, ms=4)
        plt.xticks(ticks=[i* offset_inter + ((nb_groups-1)/2 * offset_intra) for i in range(len(L_lengths))],
                   labels=L_lengths)
        plt.xlabel("Theoretical homopolymer length (reference genome)")
        plt.ylabel("Errorless sequenced homopolymer ratio (reads)")

        # Reorder legend labels
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ordered_label_list = [elt for elt in L_ordered_species if elt in by_label]
        ordered_label_values = [by_label[k] for k in ordered_label_list]
        plt.legend(ordered_label_values, ordered_label_list, title="Species:")


    plt.ylim(bottom=-1)
    plt.savefig(OUTPUT_PLOT + "percentage_homopolymer_correctly_sequenced_by_length.png")
    plt.close()

    # Save raw results in .txt file
    RAW_OUTPUT_FILE = open(OUTPUT_RAW + "percentage_homopolymer_correctly_sequenced_by_length.txt", "w")
    for group_name in dict_prct_correct:
        RAW_OUTPUT_FILE.write(f"{group_name} ({','.join([elt for elt in L_species])})\n")
        RAW_OUTPUT_FILE.write("\t".join(["Homopolymer length",
                                         "List of ratio of correctly sequenced homopolymer"]) + "\n")
        for homopolymer_length in dict_prct_correct[group_name]:
            list_values = []
            for result in dict_prct_correct[group_name][homopolymer_length]:
                if result[1] != 0:
                    list_values += [round(result[0] / result[1] * 100, 2)]
            RAW_OUTPUT_FILE.write(str(homopolymer_length) + "\t" +
                                  ", ".join(map(str, list_values))
                                  + "\n")
        RAW_OUTPUT_FILE.write("\n\n")
    RAW_OUTPUT_FILE.close()


    # 3-bis/ Modified version with details for A, C, G and T homopolymers
    # --- In case groups were defined by user: boxplot ---
    if os.path.exists(FILE_SPECIES_GROUP):
        position = 0
        nb_bases = len(L_bases)
        offset_intra = 0.25 # margin between bases (same length)
        offset_inter = (nb_bases+2) * offset_intra  # margin between length
        L_patches = []
        for i_grp, group_name in enumerate(dict_prct_correct_detail_bases):
            if dict_prct_correct_detail_bases[group_name][L_bases[0]][MIN_HOMOPOLYMER_LENGTH] == []:
                continue
            for i_base, base in enumerate(dict_prct_correct_detail_bases[group_name]):
                position = i_base * offset_intra
                color = dict_bases_colors[base]
                L_positions = [] # Position of plot
                L_prct_correct = [] # Values to plot
                for length in dict_prct_correct_detail_bases[group_name][base]:
                    tmp_list = []
                    for i, _ in enumerate(dict_prct_correct_detail_bases[group_name][base][length]):
                        if dict_prct_correct_detail_bases[group_name][base][length][i][1] != 0:
                            tmp_list += [(dict_prct_correct_detail_bases[group_name][base][length][i][0] /
                                          dict_prct_correct_detail_bases[group_name][base][length][i][1] * 100)]
                    if tmp_list == []:
                        break
                    L_prct_correct += [tmp_list]
                    L_positions += [position]
                    position += offset_inter
                box = plt.boxplot(L_prct_correct, positions=L_positions,
                                  widths=BOXPLOT_WIDTH, patch_artist=True)
                for item in ['medians']:
                    plt.setp(box[item], color="black")
                for patch in box['boxes']:
                    patch.set(facecolor=color)
                patch = mpatches.Patch(color=color, label=base)
                L_patches += [patch]
            plt.xticks(ticks=[i*offset_inter + ((nb_bases-1)/2 * offset_intra) for i in range(len(L_lengths))],
                       labels=L_lengths)
            plt.xlabel("Theoretical homopolymer length (reference genome)")
            plt.ylabel("Errorless sequenced homopolymer ratio (reads)")
            plt.legend(title="Bases:", handles=L_patches)
            plt.savefig(OUTPUT_PLOT + f"percentage_homopolymer_correctly_sequenced_by_length_base_detail_{group_name}.png")
            plt.close()

    # --- In case no groups were defined by user: plot ---
    else:
        position = 0
        nb_bases = len(L_bases)
        offset_intra = 0.401 # margin between bases (same length)
        offset_inter = (nb_bases+2) * offset_intra  # margin between length
        for i_grp, group_name in enumerate(L_ordered_species):
            if dict_prct_correct_detail_bases[group_name][L_bases[0]][MIN_HOMOPOLYMER_LENGTH] == []:
                continue
            for i_base, base in enumerate(dict_prct_correct_detail_bases[group_name]):
                position = i_base * offset_intra
                color = dict_bases_colors[base]
                L_positions = [] # Position of plot
                L_prct_correct = [] # Values to plot
                for length in dict_prct_correct_detail_bases[group_name][base]:
                    tmp_list = []
                    for i, _ in enumerate(dict_prct_correct_detail_bases[group_name][base][length]):
                        if dict_prct_correct_detail_bases[group_name][base][length][i][1] != 0:
                            tmp_list += [(dict_prct_correct_detail_bases[group_name][base][length][i][0] /
                                          dict_prct_correct_detail_bases[group_name][base][length][i][1] * 100)]
                    if tmp_list == []:
                        break
                    L_prct_correct += [tmp_list]
                    L_positions += [position]
                    position += offset_inter
                plt.plot(L_positions, L_prct_correct, "o", color=color, label=base, ms=4)
            plt.xticks(ticks=[i*offset_inter + ((nb_bases-1)/2 * offset_intra) for i in range(len(L_lengths))],
                       labels=L_lengths)
            plt.xlabel("Theoretical homopolymer length (reference genome)")
            plt.ylabel("Errorless sequenced homopolymer ratio (reads)")
            plt.legend(title="Bases:")
            plt.savefig(OUTPUT_PLOT + f"percentage_homopolymer_correctly_sequenced_by_length_base_detail_{group_name}.png")
            plt.close()
    
    # Save raw results in .txt file
    RAW_OUTPUT_FILE = open(OUTPUT_RAW + "percentage_homopolymer_correctly_sequenced_by_length_detail.txt", "w")
    for group_name in dict_prct_correct_detail_bases:
        RAW_OUTPUT_FILE.write(f"{group_name} ({','.join([elt for elt in L_species])})\n")
        RAW_OUTPUT_FILE.write("\t".join(["Homopolymer length",
                                         "List of ratio of correctly sequenced homopolymer"]) + "\n")
        for base in sorted(dict_prct_correct_detail_bases[group_name]):
            RAW_OUTPUT_FILE.write(base + "\n")
            for homopolymer_length in dict_prct_correct_detail_bases[group_name][base]:
                list_values = []
                for result in dict_prct_correct_detail_bases[group_name][base][homopolymer_length]:
                    if result[1] != 0:
                        list_values += [round(result[0] / result[1] * 100, 2)]
                RAW_OUTPUT_FILE.write(str(homopolymer_length) + "\t" +
                                      ", ".join(map(str, list_values))
                                      + "\n")
            RAW_OUTPUT_FILE.write("\n")
        RAW_OUTPUT_FILE.write("\n\n")
    RAW_OUTPUT_FILE.close()


    # 4/ --- Plot dictionary storing homopolymer length differences ---
    #     between genomic expected one, and the actual sequenced one
    
    # Set positions
    offset = 0.3
    dict_positions = {grp:[] for grp in range(len(L_groups))}
    pos = 0
    for value in range(MIN_HOMOPOLYMER_LENGTH, MAX_HOMOPOLYMER_LENGTH+1):
        for i_grp in range(len(L_groups)):
            dict_positions[i_grp] += [round(pos,1)]
            pos += offset
        pos += offset


    fig, ax = plt.subplots()
    L_patches = []
    for i_grp, group_name in enumerate(L_ordered_species):
        color = dict_species_group_color[group_name]
        if dict_diff_length_sequenced[group_name][MIN_HOMOPOLYMER_LENGTH] == {}:
            position += offset
            continue
        L_to_plot = []
        for homopolymer_genome_length in dict_diff_length_sequenced[group_name]:
            L_to_plot += [get_statistics(dict_diff_length_sequenced[group_name][homopolymer_genome_length])]

        boxprops = dict(linewidth=1)
        box = ax.bxp(L_to_plot, showfliers=False,
                     positions=dict_positions[i_grp],
                     patch_artist=True, widths=BOXPLOT_WIDTH, boxprops=boxprops)
        for item in ['boxes', 'fliers']:
            plt.setp(box[item], color=color)
        for patch in box['boxes']:
            patch.set(facecolor=color)
        for item in ['medians', 'whiskers', 'caps']:
            plt.setp(box[item], color="black")
        position += offset
        patch = mpatches.Patch(color=color, label=group_name)
        L_patches += [patch]
    # Add expected distribution x = y
    L_x, L_y = [], []
    posx = 0 
    for i in range(MIN_HOMOPOLYMER_LENGTH, MAX_HOMOPOLYMER_LENGTH + 1):
        for N in range(len(L_groups)):
            L_x += [posx]
            posx += offset      
            L_y += [i]
        posx += offset
            

        
        
    
    plt.plot(L_x, L_y, "--", color="black")
    # Details of the plot
    plt.yscale("log")
    L_tick_pos = []
    pos = 0 
    for i in range(MIN_HOMOPOLYMER_LENGTH, MAX_HOMOPOLYMER_LENGTH + 1):
        L_tick_pos += [pos + (len(L_groups)-1)/2*offset]
        pos += offset * (len(L_groups)+1)
    plt.xticks(ticks=L_tick_pos,
               labels=[length for length in dict_diff_length_sequenced[group_name]])
    plt.xlabel("Theoretical homopolymer lengths (reference genome)")
    plt.ylabel("Sequenced hompolymer lengths (reads)")
    if os.path.exists(FILE_SPECIES_GROUP):
        title = "Groups:"
    else:
        title = "Species:"
    plt.legend(handles=L_patches, title=title)
    plt.savefig(OUTPUT_PLOT + "difference_expected_sequenced_homopolymer_length.png")
    plt.close()

    # Save raw results in .txt file
    RAW_OUTPUT_FILE = open(OUTPUT_RAW + "difference_expected_sequenced_homopolymer_length.txt", "w")
    for group_name in dict_diff_length_sequenced:
        RAW_OUTPUT_FILE.write(f"{group_name}\n")
        RAW_OUTPUT_FILE.write("\t".join(["Homopolymer length",
                                         "Minimum", "Q1", "Median", "Q3", "Maximum"]) + "\n")
        if dict_diff_length_sequenced[group_name][MIN_HOMOPOLYMER_LENGTH] == {}:
            RAW_OUTPUT_FILE.write("\n\n")
            continue
        L_to_plot = []
        for homopolymer_length_genome in dict_diff_length_sequenced[group_name]:
            L_to_plot += [get_statistics(dict_diff_length_sequenced[group_name][homopolymer_length_genome])]

        for homopolymer_length_genome in dict_diff_length_sequenced[group_name]:
            result = get_statistics(dict_diff_length_sequenced[group_name][homopolymer_length_genome])
            if result["whislo"] == 0.9: # artefact for plotting log-scaled
                result["whislo"] = 0
            RAW_OUTPUT_FILE.write(str(homopolymer_length) + "\t" +
                                  str(result["whislo"]) + "\t" +
                                  str(result["q1"]) + "\t" +
                                  str(result["med"]) + "\t" +
                                  str(result["q3"]) + "\t" +
                                  str(result["whishi"]) + "\n")
        RAW_OUTPUT_FILE.write("\n\n")
    RAW_OUTPUT_FILE.close()




# --------------------------------------------------------------------------------------------------
# Main

if __name__ == "__main__":

    ## Parses arguments
    NUMBER_EXPECTED_ARGUMENTS = 9
    if len(sys.argv) != NUMBER_EXPECTED_ARGUMENTS + 1:
        print(f"ERROR: Wrong number of arguments: {NUMBER_EXPECTED_ARGUMENTS} expected but {len(sys.argv)-1} given.")
        sys.exit(2)
    ALN_EXPL_DIRNAME, REFERENCE_GENOME_DIRNAME, OUTPUT_RAW, OUTPUT_PLOT, FILE_SPECIES_GC, FILE_SPECIES_GROUP, MIN_HOMOPOLYMER_LENGTH, MAX_HOMOPOLYMER_LENGTH, NB_MAX_ALN = sys.argv[1:]
    MIN_HOMOPOLYMER_LENGTH = int(MIN_HOMOPOLYMER_LENGTH)
    MAX_HOMOPOLYMER_LENGTH = int(MAX_HOMOPOLYMER_LENGTH)
    NB_MAX_ALN = int(NB_MAX_ALN)

    if NB_MAX_ALN == -1:
        NB_MAX_ALN = float("inf")

    # --- Parameters ---

    # Colors for the graphs
    COLOR_MISMATCH = "#edd192"
    COLOR_INSERTION = "#daa4ad"
    COLOR_DELETION = "#a4beda"

    dict_color_base = {"A": '#6DA6CC', "C": '#A1DB73', "G": '#F6B255', "T": '#DC212C'}


    # Length of homopolymers
    L_lengths = list(np.arange(MIN_HOMOPOLYMER_LENGTH, MAX_HOMOPOLYMER_LENGTH)) \
                + [f"{MAX_HOMOPOLYMER_LENGTH}+"]

    # Width of boxplot
    BOXPLOT_WIDTH = 0.18

    # list of genomic bases and associated colors
    L_bases = ["A", "C", "G", "T"]
    dict_bases_colors = {"A": '#64a9ce', "C":'#a2d973', "G":"#fcb14f", "T":"#e41b1e"}

    # list of error types
    L_errors = ["Mismatch", "Insertion", "Deletion", "Total"]

    # Group species:
    #   - by user defined group (if exists)
    if os.path.exists(FILE_SPECIES_GROUP):
        dict_species_group, dict_species_group_color = get_groups(FILE_SPECIES_GROUP)
    #   - else build artificial groups: one per species
    else:
        dict_species_group, dict_species_group_color, L_ordered_species = build_group_per_species(FILE_SPECIES_GC)
    L_groups = sorted(set(dict_species_group.values()))

    # Pattern for regular expression finding of homopolymers
    pattern_homopolymer = "A{2,}|C{2,}|G{2,}|T{2,}"


    # --- Initiate dictionaries that will store results ---

    # Store genomic distribution of homopolymers, for each group
    dict_genomic_distrib = {}
    L_genomic_length = list(np.arange(MIN_HOMOPOLYMER_LENGTH, MAX_HOMOPOLYMER_LENGTH)) \
                       + [f"{MAX_HOMOPOLYMER_LENGTH}+"]
    for group_name in L_groups:
        dict_genomic_distrib[group_name] = {}
        for length in L_genomic_length:
            dict_genomic_distrib[group_name][length] = {}

    # Store number of all errors, and those linked to homopolymers
    dict_quantify_errors = {}
    for group_name in L_groups:
        dict_quantify_errors[group_name] = {"All_errors": {}, "Homopolymer_errors": {}}
        for error_type in L_errors:
            for category in dict_quantify_errors[group_name]:
                dict_quantify_errors[group_name][category][error_type] = 0


    # Store detailed errors in homopolymeric regions, depending on the error type,
    #   for each species category
    dict_error_length = initiate_dictionary_error_length()

    # Initiate dictionary that will store number of correctly sequenced homopolymers over
    #  total number of homopolymers, depending on length of genomic homopolymer
    # And a second one (same but with details for all bases A,C,G and T)
    dict_prct_correct, dict_prct_correct_detail_bases = initiate_prct_correct_dict()


    # Initiate dictionary that will store, for each genomic homopolymer's length
    #  the length of sequenced homopolymer
    dict_diff_length_sequenced = initiate_diff_len_dict()



    # --- ---
    L_species = []
    for aln_filename in os.listdir(ALN_EXPL_DIRNAME):

        aln_path = ALN_EXPL_DIRNAME + aln_filename
        if not os.path.isfile(aln_path):
            continue


        # Get total number of alignments in this file
        #     divided by 3 as alignments are represented as triplets of lines
        NB_TOT_ALN = get_total_number_of_lines(aln_path) / 3
        STARTING_TIME = time.time()

        species_name = aln_filename.split(".txt")[0]
        aln_file = open(aln_path, "r")
        group_name = dict_species_group[get_short_name(species_name)]
        L_species += [species_name]
        print("  ", aln_filename)


        # initialize dictionary that store ratio of correctly sequenced homopolymers
        for length in L_lengths:
            dict_prct_correct[group_name][length] += [[0, 0]]
            for base in L_bases:
                dict_prct_correct_detail_bases[group_name][base][length] += [[0, 0]]


        # 1/ Get homopolymer distribution from reference genome
        reference_genome_filename = REFERENCE_GENOME_DIRNAME + species_name + ".fasta"
        get_genomic_homopolymer_distribution(reference_genome_filename, group_name)


        # Get alignment from aln_file
        nb_aln_done = 0
        progressing = 0
        while True:
            header_aln = aln_file.readline().replace("\n", "")
            if not header_aln:
                break
            genome_aln = aln_file.readline().replace("\n", "")
            read_aln = aln_file.readline().replace("\n", "")

            nb_aln_done += 1
            tmp_progressing = int(nb_aln_done / NB_TOT_ALN * 100)
            if tmp_progressing > progressing and tmp_progressing % 5 == 0:
                progressing = tmp_progressing
                display_progressing_bar(progressing, time.time())
                compute_results()

            if nb_aln_done >= NB_MAX_ALN: # Limit number of alignments to speed up computation
                compute_results()
                break

            # 2/ Compute global errors for the given alignment
            update_errors(genome_aln, read_aln, group_name, "All_errors")

            # Compute sequencing errors associated with genomic homopolymers
            start_all, end_all = -1, -1
            # find homopolymers in reference genome...
            for res in re.finditer(pattern_homopolymer, genome_aln):
                start, end = res.span()
                if start in range(start_all, end_all): # this homopolymer already has been analysed
                    continue
                base_homopolymer = res.group()[0]
                start_genome, end_genome = try_extend_genomic_homopolymer(start, end, genome_aln)
                homopolymer_genome = genome_aln[start_genome: end_genome]
                homopolymer_genome_length = len(homopolymer_genome.replace("-", ""))
                if homopolymer_genome_length >= MAX_HOMOPOLYMER_LENGTH:
                    homopolymer_genome_length = f"{MAX_HOMOPOLYMER_LENGTH}+"

                # ...get associated sequenced bases
                start_read, end_read = try_extend_read_homopolymer(start, end, read_aln)
                start_all = min(start_genome, start_read)
                end_all = max(end_genome, end_read)

                hompolymer_genome_part = genome_aln[start_all: end_all]
                homopolymer_read_part = read_aln[start_all: end_all]

                # 3/ Quantification of homopolymer errors
                #    i.e. number of mismatches and indel that are linked with homopolymeric
                #    regions in the genome
                update_errors(hompolymer_genome_part, homopolymer_read_part,
                              group_name, "Homopolymer_errors")

                # 4/ Computes error rates (mismatches and indel) depending on homopolymer length
                update_errors_homopolymer_length(hompolymer_genome_part, homopolymer_read_part,
                                                 homopolymer_genome_length)

                # 5/ Update counters of correctly sequenced homopolymers, and total number
                update_prct_correct_homopolymer(hompolymer_genome_part, homopolymer_read_part,
                                                group_name,
                                                homopolymer_genome_length)

                # 5-bis/ Same for detailed dictionary
                update_prct_correct_homopolymer_detail(hompolymer_genome_part, homopolymer_read_part,
                                                       group_name, homopolymer_genome_length,
                                                       base_homopolymer)

                # 6/ Update dictionary storing homopolymer length differences
                #     between genomic expected one, and the actual sequenced one
                update_diff_len_dict(hompolymer_genome_part, homopolymer_read_part,
                                     group_name, homopolymer_genome_length, base_homopolymer)


        aln_file.close()
        sys.stdout.write("\n")


    compute_results()

