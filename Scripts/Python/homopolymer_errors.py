#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Analyses errors in homopolymeric regions, and computes (and plots) several information:
    - quantifying homopolymers in reference genomes
    - quantifying errors related to homopolymeric regions (compared to all errors)
    - error rates (mismatch, insertion, deletion) depending on homopolymer length
    - differences between expected and sequenced homopolymer lengths
All plots are made here by distinguishing low/high gc bacteria and human datasets
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
          'figure.figsize': (14, 9),
          'axes.labelsize': 16,
          'axes.titlesize': 16,
          'xtick.labelsize': 16,
          'ytick.labelsize': 16}
plt.rcParams.update(PARAMS)

# --------------------------------------------------------------------------------------------------
# Functions

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
        deletion) and their length, for each species category
        Mismatches and indels are not stored explicitely, but in the following order:
            mismatch (position 0), insertion (position 1) and deletion (position 2)
    """
    dictionary = {}
    error_length_list = np.arange(2, 8)
    for species_cat in ["bacteria", "human"]:
        dictionary[species_cat] = {}
        for error_length in error_length_list:
            dictionary[species_cat][error_length] = [0, 0, 0] # mismatch, insertion, deletion
    return dictionary

def initiate_prct_correct_dict() -> dict:
    """Initiates a dictionary which will, for each species category
    and each genomic homopolymer length, store number of correctly sequenced homopolymer,
    over total number of homopolymers.
    """
    dictionary = {}
    for species_cat in list_species_categories:
        dictionary[species_cat] = {}
        for error_length in np.arange(2, 10):
            dictionary[species_cat][error_length] = [] # list of prct for each species
    return dictionary

def initiate_prct_correct_dict_detail() -> dict:
    """Same as previous one, but with details for A, C, G and T homopolymers
    """
    dictionary = {}
    for species_cat in list_species_categories:
        dictionary[species_cat] = {}
        for base in ["A", "C", "G", "T"]:
            dictionary[species_cat][base] = {}
            for error_length in np.arange(2, 10):
                dictionary[species_cat][base][error_length] = []
    return dictionary



def initiate_diff_len_dict() -> dict:
    """Add a 'sub-dictionary' for the current species, which will, for each species category
    and for each genomic homopolymer length,
    store the sequenced homopolymer length (i.e. in the read)
    """
    dictionary = {}
    for species_cat in list_species_categories:
        dictionary[species_cat] = {}
        for error_length in np.arange(2, 10):
            dictionary[species_cat][error_length] = {}
    return dictionary


def get_genomic_homopolymer_distribution():
    """Compute homopolymer distribution in the current reference genome, and
    updates the associated dictionary genomic_distribution_homopolymer_dict
    """
    with open(reference_genome_filename, "r") as reference_genome_file:
        genome = ""
        reference_genome_file.readline() # skip first line = header
        while True:
            line = reference_genome_file.readline().rstrip()
            if not line:
                break
            if line[0] == ">":
                for res in re.finditer(pattern_homopolymer, genome):
                    start, end = res.span()
                    homopolymer_length = end - start
                    if homopolymer_length >= 10:
                        homopolymer_length = "10+"
                    base = res.group()[0]
                    if homopolymer_length not in genomic_distribution_homopolymer_dict[species_category]:
                        genomic_distribution_homopolymer_dict[species_category][homopolymer_length] = {}
                    if base not in genomic_distribution_homopolymer_dict[species_category][homopolymer_length]:
                        genomic_distribution_homopolymer_dict[species_category][homopolymer_length][base] = 0
                    genomic_distribution_homopolymer_dict[species_category][homopolymer_length][base] += 1
                genome = ""
                continue
            genome += line.upper()


    for res in re.finditer(pattern_homopolymer, genome):
        start, end = res.span()
        homopolymer_length = end - start
        if homopolymer_length >= 10:
            homopolymer_length = "10+"
        base = res.group()[0]
        if homopolymer_length not in genomic_distribution_homopolymer_dict[species_category]:
            genomic_distribution_homopolymer_dict[species_category][homopolymer_length] = {}
        if base not in genomic_distribution_homopolymer_dict[species_category][homopolymer_length]:
            genomic_distribution_homopolymer_dict[species_category][homopolymer_length][base] = 0
        genomic_distribution_homopolymer_dict[species_category][homopolymer_length][base] += 1

    # Plot
    _, ax = plt.subplots()
    dict_base_colors_light = {"A": '#64a9ce', "C":'#a2d973', "G":"#fcb14f", "T":"#e41b1e"}
    dict_base_colors_dark = {"A": '#2B6788', "C":'#62A02C', "G":"#C97303", "T":"#720D0F"}
    bar_width = 0.16

    #  for human
    dict_to_plot_human = {}
    for length in genomic_distribution_homopolymer_dict["human"]:
        pos_length = list_genomic_homopolymer_lengths.index(length)
        for base in genomic_distribution_homopolymer_dict["human"][length]:
            if base not in dict_to_plot_human:
                dict_to_plot_human[base] = [0] * len(list_genomic_homopolymer_lengths)
            mean_occurrences = genomic_distribution_homopolymer_dict["human"][length][base]
            dict_to_plot_human[base][pos_length] += mean_occurrences
    r = np.arange(len(list_genomic_homopolymer_lengths))
    for i in range(len(list_bases)):
        r = [x + bar_width for x in r]
        if i == 2:
            list_xticks = list(np.asarray(r) - bar_width / 2)
        base = list_bases[i]
        ax.bar(r, dict_to_plot_human[base], color=dict_base_colors_light[base],
               edgecolor="white", width=bar_width, label=list_bases[i])

    #  for bacteria
    dict_to_plot_bacteria = {}
    for species_cat in ["low GC", "high GC"]:
        for length in genomic_distribution_homopolymer_dict[species_cat]:
            pos_length = list_genomic_homopolymer_lengths.index(length)
            for base in genomic_distribution_homopolymer_dict[species_cat][length]:
                if base not in dict_to_plot_bacteria:
                    dict_to_plot_bacteria[base] = [0] * len(list_genomic_homopolymer_lengths)
                mean_occurrences = genomic_distribution_homopolymer_dict[species_cat][length][base]
                dict_to_plot_bacteria[base][pos_length] += mean_occurrences
    r = np.arange(len(list_genomic_homopolymer_lengths))
    for i in range(len(list_bases)):
        r = [x + bar_width for x in r]
        if i == 2:
            list_xticks = list(np.asarray(r) - bar_width / 2)
        base = list_bases[i]
        ax.bar(r, dict_to_plot_bacteria[base], color=dict_base_colors_dark[base],
               edgecolor="white", width=bar_width)

    ax.set_yscale('log')
    ax.set_xlabel("Homopolymer length", fontweight="bold")
    ax.set_ylabel("Number of occurrences", fontweight="bold")
    plt.xticks(list_xticks, list_genomic_homopolymer_lengths)

    # Add custom legend items
    human_patch = mpatches.Patch(color='#d9d9d9', label='Human')
    bacterial_patch = mpatches.Patch(color='#737373', label='Bacteria')
    first_legend = plt.legend(title="Species", handles=[human_patch, bacterial_patch],
                              loc=(0.83, 0.84))
    plt.gca().add_artist(first_legend)
    plt.ylim(0.9, 4e8)
    ax.legend(ncol=4, title="Homopolymer base:", loc=(0.53, 0.84))
    plt.savefig(OUTPUT_PLOT + "homopolymer_genomic_distribution.png")
    plt.close()

    # Save raw results in .txt file
    RAW_OUTPUT_FILE = open(OUTPUT_RAW + "homopolymer_genomic_distribution.txt", "w")
    for species_cat in list_species_categories:
        RAW_OUTPUT_FILE.write(species_cat + "\n")
        RAW_OUTPUT_FILE.write("\t".join(["Homopolymer length"] + list_bases) + "\n")
        for homopolymer_length in genomic_distribution_homopolymer_dict[species_cat]:
            list_values = list(genomic_distribution_homopolymer_dict[species_cat][homopolymer_length].values())
            RAW_OUTPUT_FILE.write(str(homopolymer_length) + "\t" +
                                  "\t".join(map(str, list_values))
                                  + "\n")



def update_global_errors():
    """ Update `quantify_all_errors_dict` (storing global error rates) with current alignment
    """
    if species_category == "human":
        quantify_all_errors_dict = quantify_all_errors_dict_human
    else:
        quantify_all_errors_dict = quantify_all_errors_dict_bact

    for i, genome_base in enumerate(genome_aln):
        read_base = read_aln[i]
        if genome_base == "N":
            continue
        if genome_base == "-":
            quantify_all_errors_dict["Insertion"] += 1
        elif read_base == "-":
            quantify_all_errors_dict["Deletion"] += 1
        elif genome_base != read_base:
            quantify_all_errors_dict["Mismatch"] += 1


def update_homopolymer_errors(g: str, r: str):
    """Computes errors in alignment of g (genome) versus r (read)
    Updates global dictionary that store these results
    """
    if species_category == "human":
        quantify_homopolymer_errors_dict = quantify_homopolymer_errors_dict_human
    else:
        quantify_homopolymer_errors_dict = quantify_homopolymer_errors_dict_bact
    for i, genome_base in enumerate(g):
        read_base = r[i]
        if genome_base == "N":
            continue
        if genome_base == "-":
            quantify_homopolymer_errors_dict["Insertion"] += 1
        elif read_base == "-":
            quantify_homopolymer_errors_dict["Deletion"] += 1
        elif genome_base != read_base:
            quantify_homopolymer_errors_dict["Mismatch"] += 1


def update_errors_homopolymer_length(g: str, r: str, homopol_len: int):
    """Computes errors in alignment of g (genome) versus r (read)
    Updates global dictionary that store these results, for length homopol_len
    Mismatches and indels are not stored explicitely, but in the following order:
        mismatch (position 0), insertion (position 1) and deletion (position 2)
    """
    if species_category == "human":
        category = "human"
    else:
        category = "bacteria"
    if homopol_len not in homopolymer_error_length_dict[category]:
        return
    for i, genome_base in enumerate(g):
        read_base = r[i]
        if genome_base == "N":
            continue
        if genome_base == "-":
            homopolymer_error_length_dict[category][homopol_len][1] += 1
        elif read_base == "-":
            homopolymer_error_length_dict[category][homopol_len][2] += 1
        elif genome_base != read_base:
            homopolymer_error_length_dict[category][homopol_len][0] += 1


def update_prct_correct_homopolymer(g: str, r: str, homopol_len: int):
    """Updates dictionary prct_correct_homopolymer_dict with counts of correctly sequenced
    homopolymers, and total number of homopolymers
    """
    if homopol_len not in prct_correct_homopolymer_dict[species_category]:
        return
    if g == r:
        prct_correct_homopolymer_dict[species_category][homopol_len][-1][0] += 1
    prct_correct_homopolymer_dict[species_category][homopol_len][-1][1] += 1

def update_prct_correct_homopolymer_detail(g: str, r: str, homopol_len: int, homopol_base: str):
    """Same as previous, but for the detailed version
    """
    if homopol_len not in prct_correct_homopolymer_dict_detail[species_category][homopol_base]:
        return
    if g == r:
        prct_correct_homopolymer_dict_detail[species_category][homopol_base][homopol_len][-1][0] += 1
    prct_correct_homopolymer_dict_detail[species_category][homopol_base][homopol_len][-1][1] += 1


def update_diff_len_dict(g: str, r: str, homopol_len_genome: int, homopol_base: str):
    """Updates dictionary that stores homopolymer length differences (genomic and sequenced)
    """
    if homopol_len_genome not in diff_homopol_length_sequenced_dict[species_category]:
        return
    homopol_len_read = len(max(re.compile(homopol_base + "*").findall(r)))
    if homopol_len_read not in diff_homopol_length_sequenced_dict[species_category][homopol_len_genome]:
        diff_homopol_length_sequenced_dict[species_category][homopol_len_genome][homopol_len_read] = 0
    diff_homopol_length_sequenced_dict[species_category][homopol_len_genome][homopol_len_read] += 1


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
    """ If sequence[start: end] start/end by '-' or homopolymer letter, extend the
    start and end positions (delimitiong genomic homopolymer)
    """
    if sequence[start_pos] in ("-", letter_homopolymer):
        while start_pos-1 >= 0 and sequence[start_pos-1] in ("-", letter_homopolymer):
            start_pos -= 1
    if sequence[end_pos-1] in ("-", letter_homopolymer):
        while end_pos < len(sequence) and sequence[end_pos] in ("-", letter_homopolymer):
            end_pos += 1
    while sequence[start_pos] == "-":
        start_pos += 1
    while sequence[end_pos-1] == "-":
        end_pos -= 1
    return start_pos, end_pos


def try_extend_read_homopolymer(start_pos, end_pos, sequence):
    """ If sequence[start: end] start/end by '-' or homopolymer letter, extend the
    start and end positions (delimitiong genomic homopolymer)
    """
    if sequence[start_pos: end_pos].replace("-", "") == "":
        return start_pos, end_pos

    if sequence[start_pos] == letter_homopolymer:
        while start_pos-1 >= 0 and sequence[start_pos-1] == letter_homopolymer:
            start_pos -= 1
    if sequence[end_pos-1] == letter_homopolymer:
        while end_pos < len(sequence) and sequence[end_pos] == letter_homopolymer:
            end_pos += 1

    return start_pos, end_pos


def compute_results():
    """Computes all results, and output as graph and raw file
    """

    # 1/ Global quantification of homopolymer errors
    #    Outputs raw .txt file + LaTeX code part to make a table
    RAW_OUTPUT_FILE = open(OUTPUT_RAW + "quantify_errors_homopolymer_vs_total.txt", "w")
    ## Raw
    categories = ["bacteria", "human"]
    for i, quantify_homopolymer_errors_dict in enumerate([quantify_homopolymer_errors_dict_bact,
                                                          quantify_homopolymer_errors_dict_human]):
        quantify_all_errors_dict = [quantify_all_errors_dict_bact, quantify_all_errors_dict_human][i]
        RAW_OUTPUT_FILE.write("---\nResults for " + categories[i] + "\n")

        error_list = ["Mismatch", "Insertion", "Deletion"]
        header = "\t".join([""] + error_list + ["Global"])
        res_homopol_err = []
        res_all_err = []
        ratio_homopol_err = []
        for error in error_list:
            nb_homopol_err = round(quantify_homopolymer_errors_dict[error] / 1e6, 2)
            nb_all_err = round(quantify_all_errors_dict[error] / 1e6, 2)
            res_homopol_err += [nb_homopol_err]
            res_all_err += [nb_all_err]
            if nb_all_err == 0.:
                ratio_homopol_err += [0]
            else:
                ratio_homopol_err += [round(nb_homopol_err / nb_all_err * 100, 2)]
        homopol_err_str = "\t".join(["H: Homopol. (x10^6)"] + \
                                     list(map(str, res_homopol_err)) + \
                                     [str(round(np.mean(res_homopol_err), 2))])
        all_err_str = "\t".join(["A: All errors (x10^6)"] + \
                                 list(map(str, res_all_err)) + \
                                 [str(round(np.mean(res_all_err), 2))])
        ratio_err = "\t".join(["Ratio H/A (%)"] + \
                               list(map(str, ratio_homopol_err)) + \
                               [str(round(np.mean(ratio_homopol_err), 2))])
        RAW_OUTPUT_FILE.write("\n".join([header, homopol_err_str, all_err_str, ratio_err]))

        ## LaTeX code
        RAW_OUTPUT_FILE.write("\n\nLaTeX code:\n\n")
        RAW_OUTPUT_FILE.write("\\begin{table}[] \n")
        RAW_OUTPUT_FILE.write("    \\centering\n")
        RAW_OUTPUT_FILE.write("    \\begin{tabular}{ccccc}\n")
        RAW_OUTPUT_FILE.write("     & Mismatches & Insertions & Deletions & Global \\\\ \n")
        RAW_OUTPUT_FILE.write("    H: Homopol. ($\\times10^{6}$) & ")
        RAW_OUTPUT_FILE.write("    " +" & ".join(list(map(str, res_homopol_err)) + \
                                         [str(round(np.mean(res_homopol_err), 2))]) + " \\\\ \n")
        RAW_OUTPUT_FILE.write("    A: All errors ($\\times10^{6}$) &")
        RAW_OUTPUT_FILE.write("    " + " & ".join(list(map(str, res_all_err)) + \
                                         [str(round(np.mean(res_all_err), 2))]) + "  \\\\ \n")
        RAW_OUTPUT_FILE.write("    Ratio H/A (\\%) &\n")
        RAW_OUTPUT_FILE.write("    " + " & ".join(list(map(str, ratio_homopol_err)) + \
                                         [str(round(np.mean(ratio_homopol_err), 2))]) + "  \\\\ \n")
        RAW_OUTPUT_FILE.write("    \\end{tabular}\n")
        RAW_OUTPUT_FILE.write("    \\caption{Counts of total and homopolymer-induced sequencing errors.}\n")
        RAW_OUTPUT_FILE.write("\\end{table}\n\n\n")

    RAW_OUTPUT_FILE.close()


    # 2/ Computes a pie of error rates (mismatches and indel) depending on homopolymer length
    # 1) save a plot
    fig = plt.figure()
    text_size = 16
    pos_x = 0
    pos_y = 0.52
    color_list = [COLOR_MISMATCH, COLOR_INSERTION, COLOR_DELETION]
    label_list = ["Mismatches", "Insertions", "Deletions"]

    #    bacteria
    for length in homopolymer_error_length_dict["bacteria"]:
        ax = fig.add_axes([pos_x, pos_y, .5, .5], aspect=1)
        wedges, _, _ = ax.pie(homopolymer_error_length_dict["bacteria"][length],
                              colors=color_list, autopct='%1.0f%%',
                              textprops=dict(color="k", size=text_size),
                              pctdistance=0.8,
                              radius=0.85,
                              startangle=90)
        for w in wedges:
            w.set_linewidth(2)
            w.set_edgecolor('k')
        title = "Length " + str(length)
        ax.set_title(title, y=0.07)
        pos_x += 0.25
        if pos_x > 0.5:
            pos_x = 0
            pos_y = 0
        fig.legend(loc=10, title="Error:", ncol=3, labels=label_list)

    #    human
    pos_x = 0
    pos_y = 0.52
    for length in homopolymer_error_length_dict["human"]:
        ax = fig.add_axes([pos_x, pos_y, .5, .5], aspect=1)
        wedges, _, _ = ax.pie(homopolymer_error_length_dict["human"][length],
                              colors=color_list, autopct='%1.0f%%',
                              textprops=dict(color="k", size=text_size),
                              pctdistance=0.7,
                              radius=0.5,
                              startangle=90)
        for w in wedges:
            w.set_linewidth(2)
            w.set_edgecolor('k')
        title = "Length " + str(length)
        ax.set_title(title, y=0.07)
        pos_x += 0.25
        if pos_x > 0.5:
            pos_x = 0
            pos_y = 0

    plt.savefig(OUTPUT_PLOT + "error_rates_depending_on_homopolymer_length.png")
    plt.close()


    # 2) save raw results in .txt file
    RAW_OUTPUT_FILE = open(OUTPUT_RAW + "error_rates_depending_on_homopolymer_length.txt", "w")
    for category in ["bacteria", "human"]:
        RAW_OUTPUT_FILE.write("---\nResults for " + category + "\n")

        #   Raw
        error_list = ["Mismatch", "Insertion", "Deletion"]
        RAW_OUTPUT_FILE.write("\t".join(["Homopolymer length"] + error_list) + "\n")
        for homopolymer_length in homopolymer_error_length_dict[category]:
            RAW_OUTPUT_FILE.write(str(homopolymer_length) + "\t" +
                                  "\t".join(map(str, homopolymer_error_length_dict[category][homopolymer_length]))
                                  + "\n")
        #   LaTeX code
        RAW_OUTPUT_FILE.write("\n\nLaTeX code:\n\n")
        RAW_OUTPUT_FILE.write("\\begin{table}[] \n")
        RAW_OUTPUT_FILE.write("    \\centering\n")
        RAW_OUTPUT_FILE.write("    \\begin{tabular}{cccc}\n")
        RAW_OUTPUT_FILE.write("    " + " & ".join(["Homopolymer length"] + error_list) + " \\\\ \n")
        for homopolymer_length in homopolymer_error_length_dict[category]:
            RAW_OUTPUT_FILE.write(" " * 4 + str(homopolymer_length) + " & " +
                                  " & ".join(map(str, homopolymer_error_length_dict[category][homopolymer_length]))
                                  + " \\\\ \n")
        RAW_OUTPUT_FILE.write("    \\end{tabular}\n")
        RAW_OUTPUT_FILE.write("    \\caption{Error counts depending on homopolymer sizes, for " +
                                             "each error type.}\n")
        RAW_OUTPUT_FILE.write("\\end{table}\n\n\n")

    RAW_OUTPUT_FILE.close()



    # 3/ Plots ratio of correctly sequenced homopolymers over and total number of homopolymers
    # Prepare data for the plot, and converts counts to ratios
    offset = -0.25
    list_color = [COLOR_LOW_GC, COLOR_HIGH_GC, COLOR_HUMAN]
    for i_category, category in enumerate(prct_correct_homopolymer_dict):
        lengths_list = []
        prct_correct_list = []
        color = list_color[i_category]
        for length in prct_correct_homopolymer_dict[category]:
            lengths_list += [length + offset]
            tmp_list = []
            for i, _ in enumerate(prct_correct_homopolymer_dict[category][length]):
                if prct_correct_homopolymer_dict[category][length][i][1] != 0:
                    tmp_list += [(prct_correct_homopolymer_dict[category][length][i][0] /
                                  prct_correct_homopolymer_dict[category][length][i][1] * 100)]
            prct_correct_list += [tmp_list]
        # Plot
        box = plt.boxplot(prct_correct_list, positions=lengths_list,
                          widths=BOXPLOT_WIDTH, patch_artist=True)
        for item in ['medians']:
            plt.setp(box[item], color="black")
        for patch in box['boxes']:
            patch.set(facecolor=color)
        offset += 0.25
    plt.xticks(ticks=[length for length in prct_correct_homopolymer_dict[category]],
               labels=[length for length in prct_correct_homopolymer_dict[category]])
    plt.xlabel("Theoretical homopolymer length (reference genome)")
    plt.ylabel("Errorless sequenced homopolymer ratio (reads)")

    # Add custom legend items
    human_patch = mpatches.Patch(color=COLOR_HUMAN, label='Human')
    bacterial_low_patch = mpatches.Patch(color=COLOR_LOW_GC, label='Bacteria low GC')
    bacterial_high_patch = mpatches.Patch(color=COLOR_HIGH_GC, label='Bacteria high GC')
    first_legend = plt.legend(title="Species:", handles=[bacterial_low_patch,
                                                        bacterial_high_patch,
                                                        human_patch])
    plt.gca().add_artist(first_legend)
    plt.xlim(1.5, 9.5)
    plt.ylim(-1, 100)
    plt.yticks(ticks=np.arange(0, 101, 10), labels=np.arange(0, 101, 10))
    plt.savefig(OUTPUT_PLOT + "percentage_homopolymer_correctly_sequenced_by_length.png")
    plt.close()

    # Save raw results in .txt file
    RAW_OUTPUT_FILE = open(OUTPUT_RAW + "percentage_homopolymer_correctly_sequenced_by_length.txt", "w")
    for category in prct_correct_homopolymer_dict:
        RAW_OUTPUT_FILE.write("---\nResults for " + category + "\n")
        RAW_OUTPUT_FILE.write("\t".join(["Homopolymer length",
                                         "List of ratio of correctly sequenced homopolymer"]) + "\n")
        for homopolymer_length in prct_correct_homopolymer_dict[category]:
            list_values = []
            for result in prct_correct_homopolymer_dict[category][homopolymer_length]:
                if result[1] != 0:
                    list_values += [round(result[0] / result[1] * 100, 2)]
            RAW_OUTPUT_FILE.write(str(homopolymer_length) + "\t" +
                                  ", ".join(map(str, list_values))
                                  + "\n")
        RAW_OUTPUT_FILE.write("\n\n")
    RAW_OUTPUT_FILE.close()


    # 3-bis/ Modified version with details for A, C, G and T homopolymers
    spacing_between_length = 1.42

    # Add custom legend items
    A_patch = mpatches.Patch(color='#6DA6CC', label='A')
    C_patch = mpatches.Patch(color='#A1DB73', label='C')
    G_patch = mpatches.Patch(color='#F6B255', label='G')
    T_patch = mpatches.Patch(color='#DC212C', label='T')

    # Prepare data for the plot, and converts counts to ratios
    list_colors = ["#6DA6CC", "#A1DB73", "#F6B255", "#DC212C"]
    # for low gc bacteria
    offset = -0.375
    category = "low GC"
    for i, base in enumerate(prct_correct_homopolymer_dict_detail[category]):
        lengths_list = []
        prct_correct_list = []
        color = list_colors[i]
        for length in prct_correct_homopolymer_dict_detail[category][base]:
            lengths_list += [length*spacing_between_length + offset]
            tmp_list = []
            for i, _ in enumerate(prct_correct_homopolymer_dict_detail[category][base][length]):
                if prct_correct_homopolymer_dict_detail[category][base][length][i][1] != 0:
                    tmp_list += [(prct_correct_homopolymer_dict_detail[category][base][length][i][0] /
                                  prct_correct_homopolymer_dict_detail[category][base][length][i][1] * 100)]
            prct_correct_list += [tmp_list]
        # Plot
        box = plt.boxplot(prct_correct_list, positions=lengths_list,
                          widths=BOXPLOT_WIDTH, patch_artist=True)
        for item in ['medians']:
            plt.setp(box[item], color="black")
        for patch in box['boxes']:
            patch.set(facecolor=color)
        offset += 0.25
    plt.xticks(ticks=[length*spacing_between_length for length in prct_correct_homopolymer_dict_detail[category][base]],
               labels=[length for length in prct_correct_homopolymer_dict_detail[category][base]])
    plt.ylim(0, 100)
    plt.yticks(ticks=np.arange(0, 101, 10), labels=np.arange(0, 101, 10))
    plt.xlim(2, 14)
    plt.xlabel("Theoretical homopolymer length (reference genome)")
    plt.ylabel("Errorless sequenced homopolymer ratio (reads)")
    first_legend = plt.legend(title="Bases:", handles=[A_patch, C_patch, G_patch, T_patch], ncol=4)
    plt.gca().add_artist(first_legend)
    plt.savefig(OUTPUT_PLOT + "percentage_homopolymer_correctly_sequenced_by_length_detail_low_GC_bacteria.png")
    plt.close()


    # for high gc bacteria
    offset = -0.375
    category = "high GC"
    for i, base in enumerate(prct_correct_homopolymer_dict_detail[category]):
        lengths_list = []
        prct_correct_list = []
        color = list_colors[i]
        for length in prct_correct_homopolymer_dict_detail[category][base]:
            lengths_list += [length*spacing_between_length + offset]
            tmp_list = []
            for i, _ in enumerate(prct_correct_homopolymer_dict_detail[category][base][length]):
                if prct_correct_homopolymer_dict_detail[category][base][length][i][1] != 0:
                    tmp_list += [(prct_correct_homopolymer_dict_detail[category][base][length][i][0] /
                                  prct_correct_homopolymer_dict_detail[category][base][length][i][1] * 100)]
            prct_correct_list += [tmp_list]
        # Plot
        box = plt.boxplot(prct_correct_list, positions=lengths_list,
                          widths=BOXPLOT_WIDTH, patch_artist=True)
        for item in ['medians']:
            plt.setp(box[item], color="black")
        for patch in box['boxes']:
            patch.set(facecolor=color)
        offset += 0.25
    plt.xticks(ticks=[length*spacing_between_length for length in prct_correct_homopolymer_dict_detail[category][base]],
               labels=[length for length in prct_correct_homopolymer_dict_detail[category][base]])
    plt.ylim(0, 100)
    plt.yticks(ticks=np.arange(0, 101, 10), labels=np.arange(0, 101, 10))
    plt.xlabel("Theoretical homopolymer length (reference genome)")
    plt.ylabel("Errorless sequenced homopolymer ratio (reads)")
    plt.xlim(2, 14)
    first_legend = plt.legend(title="Bases:", handles=[A_patch, C_patch, G_patch, T_patch], ncol=4)
    plt.gca().add_artist(first_legend)
    plt.savefig(OUTPUT_PLOT + "percentage_homopolymer_correctly_sequenced_by_length_detail_high_GC_bacteria.png")
    plt.close()


    # for human
    offset = -0.375
    category = "human"
    for i, base in enumerate(prct_correct_homopolymer_dict_detail[category]):
        lengths_list = []
        prct_correct_list = []
        color = list_colors[i]
        for length in prct_correct_homopolymer_dict_detail[category][base]:
            lengths_list += [length*spacing_between_length + offset]
            tmp_list = []
            for i, _ in enumerate(prct_correct_homopolymer_dict_detail[category][base][length]):
                if prct_correct_homopolymer_dict_detail[category][base][length][i][1] != 0:
                    tmp_list += [(prct_correct_homopolymer_dict_detail[category][base][length][i][0] /
                                  prct_correct_homopolymer_dict_detail[category][base][length][i][1] * 100)]
            prct_correct_list += [tmp_list]
        # Plot
        box = plt.boxplot(prct_correct_list, positions=lengths_list,
                          widths=BOXPLOT_WIDTH, patch_artist=True)
        for item in ['medians']:
            plt.setp(box[item], color="black")
        for patch in box['boxes']:
            patch.set(facecolor=color)
        offset += 0.25
    plt.xticks(ticks=[length*spacing_between_length for length in prct_correct_homopolymer_dict_detail[category][base]],
               labels=[length for length in prct_correct_homopolymer_dict_detail[category][base]])
    plt.ylim(0, 100)
    plt.yticks(ticks=np.arange(0, 101, 10), labels=np.arange(0, 101, 10))
    plt.xlim(2, 14)
    plt.xlabel("Theoretical homopolymer length (reference genome)")
    plt.ylabel("Errorless sequenced homopolymer ratio (reads)")
    first_legend = plt.legend(title="Bases:", handles=[A_patch, C_patch, G_patch, T_patch], ncol=4)
    plt.gca().add_artist(first_legend)
    plt.savefig(OUTPUT_PLOT + "percentage_homopolymer_correctly_sequenced_by_length_detail_human.png")
    plt.close()



    # Save raw results in .txt file
    RAW_OUTPUT_FILE = open(OUTPUT_RAW + "percentage_homopolymer_correctly_sequenced_by_length_detail.txt", "w")
    for category in prct_correct_homopolymer_dict_detail:
        RAW_OUTPUT_FILE.write("---\nResults for " + category + "\n")
        RAW_OUTPUT_FILE.write("\t".join(["Homopolymer length",
                                         "List of ratio of correctly sequenced homopolymer"]) + "\n")
        for base in sorted(prct_correct_homopolymer_dict_detail[category]):
            RAW_OUTPUT_FILE.write(base + "\n")
            for homopolymer_length in prct_correct_homopolymer_dict_detail[category][base]:
                list_values = []
                for result in prct_correct_homopolymer_dict_detail[category][base][homopolymer_length]:
                    if result[1] != 0:
                        list_values += [round(result[0] / result[1] * 100, 2)]
                RAW_OUTPUT_FILE.write(str(homopolymer_length) + "\t" +
                                      ", ".join(map(str, list_values))
                                      + "\n")
        RAW_OUTPUT_FILE.write("\n\n")
    RAW_OUTPUT_FILE.close()







    # 4/ Update dictionary storing homopolymer length differences
    #     between genomic expected one, and the actual sequenced one
    fig, ax = plt.subplots()
    offset = -0.25
    list_colors = [COLOR_LOW_GC, COLOR_HIGH_GC, COLOR_HUMAN]
    for category in diff_homopol_length_sequenced_dict:
        if diff_homopol_length_sequenced_dict[category][2] == {}:
            offset += 0.25
            continue
        L_to_plot = []
        for homopolymer_length_genome in diff_homopol_length_sequenced_dict[category]:
            L_to_plot += [get_statistics(diff_homopol_length_sequenced_dict[category][homopolymer_length_genome])]
        COLOR = list_colors[0]
        list_colors = list_colors[1:]
        boxprops = dict(linewidth=4)
        box = ax.bxp(L_to_plot, showfliers=False,
                     positions=np.arange(MIN_HOMOPOLYMER_LENGTH + offset,
                                         MAX_HOMOPOLYMER_LENGTH + offset),
                     patch_artist=True,
                     widths=BOXPLOT_WIDTH,
                     boxprops=boxprops)
        for item in ['boxes', 'fliers']:
            plt.setp(box[item], color=COLOR)
        for patch in box['boxes']:
            patch.set(facecolor=COLOR)
        for item in ['medians', 'whiskers', 'caps']:
            plt.setp(box[item], color="black")
        offset += 0.25


    # Plots expected distribution x = y
    L_x, L_y = [], []
    for i in range(MIN_HOMOPOLYMER_LENGTH, MAX_HOMOPOLYMER_LENGTH):
        L_x += [i - 0.2, i, i+0.2]
        L_y += [i, i, i]
    plt.plot(L_x, L_y, "--", color="black")
    # Details of the plot
    plt.yscale("log")
    L_ticks = [0.9, 2, 5, 10, 20, 50, 100]
    L_ticks_labels = [0, 2, 5, 10, 20, 50, 100]
    plt.xticks(ticks=[length for length in diff_homopol_length_sequenced_dict[category]],
               labels=[length for length in diff_homopol_length_sequenced_dict[category]])
    plt.yticks(np.asarray(L_ticks), labels=L_ticks_labels)
    plt.ylim(0.8, 150)
    plt.xlim(1.5, 9.5)
    plt.xlabel("Theoretical homopolymer lengths (reference genome)")
    plt.ylabel("Errorless sequenced homopolymers ratio (reads)")

    low_patch = mpatches.Patch(color=COLOR_LOW_GC, label='Low GC bacteria')
    high_patch = mpatches.Patch(color=COLOR_HIGH_GC, label='High GC bacteria')
    human_patch = mpatches.Patch(color=COLOR_HUMAN, label='Human')
    plt.legend(handles=[low_patch, high_patch, human_patch],
               title="Species:", loc=9, ncol=3)



    plt.savefig(OUTPUT_PLOT + "difference_expected_sequenced_homopolymer_length.png")
    plt.close()


    # Save raw results in .txt file
    RAW_OUTPUT_FILE = open(OUTPUT_RAW + "difference_expected_sequenced_homopolymer_length.txt", "w")
    for category in diff_homopol_length_sequenced_dict:
        RAW_OUTPUT_FILE.write("---\nResults for " + category + "\n")
        RAW_OUTPUT_FILE.write("\t".join(["Homopolymer length",
                                         "Minimum", "Q1", "Median", "Q3", "Maximum"]) + "\n")
        if diff_homopol_length_sequenced_dict[category][2] == {}:
            RAW_OUTPUT_FILE.write("\n\n")
            continue
        L_to_plot = []
        for homopolymer_length_genome in diff_homopol_length_sequenced_dict[category]:
            L_to_plot += [get_statistics(diff_homopol_length_sequenced_dict[category][homopolymer_length_genome])]

        for homopolymer_length_genome in diff_homopol_length_sequenced_dict[category]:
            result = get_statistics(diff_homopol_length_sequenced_dict[category][homopolymer_length_genome])
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
    if len(sys.argv) != 6:
        print(f"ERROR: Wrong number of arguments: 5 expected but {len(sys.argv)-1} given.")
        sys.exit(2)
    ALN_EXPL_DIRNAME = sys.argv[1]
    REFERENCE_GENOME_DIRNAME = sys.argv[2]
    OUTPUT_RAW = sys.argv[3]
    OUTPUT_PLOT = sys.argv[4]
    FILENAME_SPECIES_GC = sys.argv[5]

    if ALN_EXPL_DIRNAME[-1] != "/":
        ALN_EXPL_DIRNAME += "/"
    if REFERENCE_GENOME_DIRNAME[-1] != "/":
        REFERENCE_GENOME_DIRNAME += "/"
    if OUTPUT_PLOT[-1] != "/":
        OUTPUT_PLOT += "/"
    if OUTPUT_RAW[-1] != "/":
        OUTPUT_RAW += "/"

    # Colors for the graphs
    COLOR_MISMATCH = "#edd192"
    COLOR_INSERTION = "#daa4ad"
    COLOR_DELETION = "#a4beda"

    COLOR_LOW_GC = "#4477AA"
    COLOR_HIGH_GC = "#EE6677"
    COLOR_HUMAN = "#CCBB44"

    MIN_HOMOPOLYMER_LENGTH, MAX_HOMOPOLYMER_LENGTH = 2, 10
    # Width of boxplot
    BOXPLOT_WIDTH = 0.2

    # Initiate dictionary that will store genomic distribution of homopolymers
    genomic_distribution_homopolymer_dict = {}
    list_genomic_homopolymer_lengths = list(np.arange(2, 10)) + ["10+"]
    list_bases = ["A", "C", "G", "T"]
    list_species_categories = ["low GC", "high GC", "human"]
    for species_cat in list_species_categories:
        genomic_distribution_homopolymer_dict[species_cat] = {}
        for length in list_genomic_homopolymer_lengths:
            genomic_distribution_homopolymer_dict[species_cat][length] = {}
            for base in list_bases:
                genomic_distribution_homopolymer_dict[species_cat][length][base] = 0


    # Initiate dictionaries that will store number of all errors, and those linked to homopolymers
    #   for bacteria
    quantify_all_errors_dict_bact = {'Mismatch': 0, 'Insertion': 0, 'Deletion': 0}
    quantify_homopolymer_errors_dict_bact = {'Mismatch': 0, 'Insertion': 0, 'Deletion': 0}
    #   for human
    quantify_all_errors_dict_human = {'Mismatch': 0, 'Insertion': 0, 'Deletion': 0}
    quantify_homopolymer_errors_dict_human = {'Mismatch': 0, 'Insertion': 0, 'Deletion': 0}

    # Initiate dictionary that will store detailed errors in homopolymeric regions, depending
    #  on the error type, for each species category
    homopolymer_error_length_dict = initiate_dictionary_error_length()

    # Initiate dictionary that will store number of correctly sequenced homopolymers over
    #  total number of homopolymers, depending on length of genomic homopolymer
    prct_correct_homopolymer_dict = initiate_prct_correct_dict()

    # Same but with details for A, C, G and T
    prct_correct_homopolymer_dict_detail = initiate_prct_correct_dict_detail()

    # Initiate dictionary that will store, for each genomic homopolymer's length
    #  the length of sequenced homopolymer
    diff_homopol_length_sequenced_dict = initiate_diff_len_dict()

    # Pattern for regular expression finding of homopolymers
    pattern_homopolymer = "A{2,}|C{2,}|G{2,}|T{2,}"

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


    for aln_filename in os.listdir(ALN_EXPL_DIRNAME):

        aln_path = ALN_EXPL_DIRNAME + aln_filename
        if not os.path.isfile(aln_path):
            continue

        NB_TOT_ALN = get_total_number_of_lines(aln_path) / 3
        STARTING_TIME = time.time()

        print(aln_filename)
        species_name = aln_filename.split(".txt")[0]
        aln_file = open(aln_path, "r")

        for species_cat in list_species_categories:
            for length in prct_correct_homopolymer_dict[species_cat]:
                prct_correct_homopolymer_dict[species_cat][length] += [[0, 0]]

        for species_cat in list_species_categories:
            for base in ["A", "C", "G", "T"]:
                for length in prct_correct_homopolymer_dict_detail[species_cat][base]:
                    prct_correct_homopolymer_dict_detail[species_cat][base][length] += [[0, 0]]

        species_category = dict_species_gc_category[species_name]

        # Get homopolymer distribution from reference genome
        reference_genome_filename = REFERENCE_GENOME_DIRNAME + species_name + ".fasta"
        get_genomic_homopolymer_distribution()

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

            #if nb_aln_done >= 100000: # Limit to 100 000 alignments
            #    print("Limit of 100 000 alignments reached. Move to next species")
            #    compute_results()
            #    break

            # Compute global errors for the given alignment
            update_global_errors()

            # Compute sequencing errors associated with genomic homopolymers
            start_all, end_all = -1, -1
            for res in re.finditer(pattern_homopolymer, genome_aln):
                start, end = res.span()
                if start in range(start_all, end_all): # this homopolymer already has been analysed
                    continue
                letter_homopolymer = res.group()[0]
                start_genome, end_genome = try_extend_genomic_homopolymer(start, end, genome_aln)
                homopolymer_genome = genome_aln[start_genome: end_genome]
                homopolymer_genome_length = len(homopolymer_genome.replace("-", ""))
                if homopolymer_genome_length > 9:
                    continue
                start_read, end_read = try_extend_read_homopolymer(start, end, read_aln)
                start_all = min(start_genome, start_read)
                end_all = max(end_genome, end_read)


                # 1/ Quantification of homopolymer errors
                #    i.e. number of mismatches and indel that are linked with homopolymeric
                #    regions in the genome
                update_homopolymer_errors(genome_aln[start_all: end_all],
                                          read_aln[start_all: end_all])


                # 2/ Computes error rates (mismatches and indel) depending on homopolymer length
                update_errors_homopolymer_length(genome_aln[start_all: end_all],
                                                 read_aln[start_all: end_all],
                                                 homopolymer_genome_length)

                # 3/ Update counters of correctly sequenced homopolymers, and total number
                update_prct_correct_homopolymer(genome_aln[start_all: end_all],
                                                read_aln[start_all: end_all],
                                                homopolymer_genome_length)


                # 3-bis/ Same for detailed dictionary
                update_prct_correct_homopolymer_detail(genome_aln[start_all: end_all],
                                                read_aln[start_all: end_all],
                                                homopolymer_genome_length,
                                                letter_homopolymer)


                # 4/ Update dictionary storing homopolymer length differences
                #     between genomic expected one, and the actual sequenced one
                update_diff_len_dict(genome_aln[start_genome: end_genome],
                                     read_aln[start_read: end_read],
                                     homopolymer_genome_length,
                                     letter_homopolymer)



        aln_file.close()
        sys.stdout.write("\n")


    compute_results()
