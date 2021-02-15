#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Analyses errors in heteropolymeric regions, and computes (and plots) several information:
    - quantifying heteropolymers in reference genomes
    - differences between expected and sequenced heteropolymer lengths
    - error rate and abundance in reads depending on heteropolymer type
"""

# --------------------------------------------------------------------------------------------------
# Packages

import sys
import os
import re
import time
import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D
# Parameters for output plots
PARAMS = {'legend.fontsize': 20,
          'legend.title_fontsize': 20,
          'legend.labelspacing': 0.1,
          'legend.borderpad': 0.3,
          'legend.columnspacing': 1,
          'legend.handletextpad': 0.5,
          'legend.handlelength': 0.8,
          'figure.figsize': (14, 9),
          'axes.labelsize': 20,
          'axes.titlesize': 22,
          'xtick.labelsize': 20,
          'ytick.labelsize': 20}
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


def initiate_diff_len_dict() -> dict:
    """Add a 'sub-dictionary' for the current species, which will, for each genomic
    heteropolymer length, store the sequenced heteropolymer length (i.e. in the read)
    """
    dictionary = {}
    for gc_category in set(dict_species_gc_category.values()):
        dictionary[gc_category] = {}
        for error_length in range(4, 12, 2):
            dictionary[gc_category][error_length] = {}
    return dictionary


def generate_heteropolymer_pattern():
    """Generates a regular expression to find heteropolymers
    i.e. repetition of groups of 2 nucleotides
    """
    list_bases = ["A", "C", "G", "T"]
    pattern = ""
    for b1 in list_bases:
        for b2 in list_bases:
            if b1 == b2:
                continue
            pattern += f"(?P<{b1+b2}>({b1+b2}){{2,}})|"
    pattern = pattern[:-1]
    return pattern

def get_genomic_heteropolymer_distribution():
    """Compute heteropolymer distribution in the current reference genome, and
    updates the associated dictionary genomic_distribution_heteropolymer_dict
    Also returns length of genome
    """
    with open(reference_genome_filename, "r") as reference_genome_file:
        genome = ""
        while True:
            line = reference_genome_file.readline().replace("\n", "")
            if not line:
                break
            if line[0] == ">":
                continue
            genome += line
    reference_genome_length = len(genome)

    for res in re.finditer(pattern_heteropolymer, genome):
        start, end = res.span()
        heteropolymer_length = end - start
        if heteropolymer_length >= 12:
            heteropolymer_length = "12+"
        base = res.group()[:2]
        # use canonical form:
        if base[::-1] < base:
            base = base[::-1]
        genomic_distribution_heteropolymer_dict[heteropolymer_length][base] += 1

    # Prepare data for plot
    bar_width = 0.1
    dict_to_plot = {}
    length_color_dict = {4: "#1965b0", 6: "#7bafde", 8: "#cae0ab", 10: "#f1932d", "12+": "#a5170e"}
    for length in genomic_distribution_heteropolymer_dict:
        dict_to_plot[length] = []
        for category in list_heteropolymer_category:
            occurrence = round(genomic_distribution_heteropolymer_dict[length][category], 2)
            dict_to_plot[length] += [occurrence]

    # Plot results:
    L_labels = ['AC-CA', 'AG-GA', 'AT-TA', 'CG-GC', 'CT-TC', 'GT-TG']
    fig, ax = plt.subplots(constrained_layout=True)
    r = np.arange(len(list_heteropolymer_category))
    for length in list_genomic_heteropolymer_lengths:
        r = [x + bar_width for x in r]
        ax.bar(r, dict_to_plot[length],
               color=length_color_dict[length], width=bar_width,
               edgecolor="white", label=length)
    ax.set_yscale('log')
    ax.set_xlabel("Heteropolymer type", fontweight="bold")
    ax.set_ylabel("Occurrences", fontweight="bold")
    plt.xticks([r + 3*bar_width for r in range(len(genomic_distribution_heteropolymer_dict[list_genomic_heteropolymer_lengths[0]]))],
               L_labels)
    plt.ylim(1e2, 5*1e8)
    ax.legend(ncol=5, title="Heteropolymer length", loc=1)
    plt.savefig(OUTPUT_PLOT + "heteropolymer_genomic_distribution.png")
    plt.close()

    # Save raw results in .txt file
    RAW_OUTPUT_FILE = open(OUTPUT_RAW + "heteropolymer_genomic_distribution.txt", "w")
    RAW_OUTPUT_FILE.write("\t".join(["Heteropolymer length"] + L_labels) + "\n")
    for heteropolymer_length in genomic_distribution_heteropolymer_dict:
        RAW_OUTPUT_FILE.write(str(heteropolymer_length))
        for category in list_heteropolymer_category:
            occurrence = genomic_distribution_heteropolymer_dict[heteropolymer_length][category]
            RAW_OUTPUT_FILE.write("\t" + str(occurrence))
        RAW_OUTPUT_FILE.write("\n")
    RAW_OUTPUT_FILE.close()

    return reference_genome_length



def update_diff_len_dict(g: str, r: str, heteropol_len_genome: int, heteropol_category: str):
    """Updates dictionary that stores heteropolymer length differences (genomic and sequenced)
    """
    g = g.replace("-", "")
    r = r.replace("-", "")
    if heteropol_len_genome not in diff_heteropol_length_sequenced_dict[gc_category]:
        return
    if heteropol_category not in r:
        heteropol_len_read = 0
    else:
        potential_heteropol_len_read = []
        pattern = "(" + heteropol_category + "){0,}"
        for match in re.finditer(pattern, r):
            match_length = match.span()[1] - match.span()[0]
            potential_heteropol_len_read += [match_length]
        heteropol_len_read = max(potential_heteropol_len_read)
    if heteropol_len_read not in diff_heteropol_length_sequenced_dict[gc_category][heteropol_len_genome]:
        diff_heteropol_length_sequenced_dict[gc_category][heteropol_len_genome][heteropol_len_read] = 0
    diff_heteropol_length_sequenced_dict[gc_category][heteropol_len_genome][heteropol_len_read] += 1


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

    # Get N = total number of "observation"
    N = 0
    for heteropol_read in dictionary:
        N += dictionary[heteropol_read]

    # Minimum, maximum, median and quartiles:
    minimum = sorted(list(dictionary))[0]
    maximum = sorted(list(dictionary))[-1]
    N_q1 = N / 4
    N_med = N / 2
    N_q3 = N * 3/4
    q1, median, q3 = None, None, None
    n0 = 0
    n1 = n0
    for heteropol_read in sorted(dictionary):
        n0 = n1
        n1 = n0 + dictionary[heteropol_read]
        # First quartile
        if q1 is None:
            if n0 < N_q1 < n1:
                q1 = heteropol_read
            elif N_q1 < n1:
                if n0 == N_q1:
                    q1 = heteropol_read
        # Median
        if median is None:
            if n0 < N_med < n1:
                median = heteropol_read
            elif N_med < n1:
                if n0 == N_med:
                    median = heteropol_read
        # Third quartile
        if q3 is None:
            if n0 < N_q3 < n1:
                q3 = heteropol_read
            elif N_q3 < n1:
                if n0 == N_q3:
                    q3 = heteropol_read
    # Turns 0 values to 1 for log-scaled output graph
    if minimum == 0:
        minimum = 0.9
    return({'med': median, 'q1': q1, 'q3': q3, 'whislo': minimum, 'whishi': maximum})


def try_extend_genomic_heteropolymer(start_pos, end_pos, sequence):
    """ If sequence[start: end] start/end by '-' or heteropolymer category, extend the
    start and end positions (delimiting genomic heteropolymer)
    """
    # Try extend left
    new_start_pos = start_pos
    while True:
        while new_start_pos - 1 >= 0 and sequence[new_start_pos - 1] == "-":
            new_start_pos -= 1
        if new_start_pos - 1 >= 0 and sequence[new_start_pos - 1] == category_heteropolymer[-1]:
            new_start_pos -= 1
        else:
            break
        while new_start_pos - 1 >= 0 and sequence[new_start_pos - 1] == "-":
            new_start_pos -= 1
        if new_start_pos - 1 >= 0 and sequence[new_start_pos - 1] == category_heteropolymer[0]:
            new_start_pos -= 1
        else:
            break
        start_pos = new_start_pos

    # Try extend right
    new_end_pos = end_pos
    while True:
        while new_end_pos + 1 < len(sequence) and sequence[new_end_pos] == "-":
            new_end_pos += 1
        if new_end_pos + 1 < len(sequence) and sequence[new_end_pos] == category_heteropolymer[0]:
            new_end_pos += 1
        else:
            break
        while new_end_pos + 1 < len(sequence) and sequence[new_end_pos] == "-":
            new_end_pos += 1
        if new_end_pos + 1 < len(sequence) and sequence[new_end_pos] == category_heteropolymer[1]:
            new_end_pos += 1
        else:
            break
        end_pos = new_end_pos

    return start_pos, end_pos



def try_extend_read_heteropolymer(start_pos, end_pos, sequence):
    """ If sequence[start: end] start/end by '-' or heteropolymer letter, extend the
    start and end positions (delimiting sequenced heteropolymer)
    Same function as try_extend_genomic_heteropolymer but here apply constraints on maximum
      gap length
    """
    if sequence[start_pos: end_pos].replace("-", "") == "":
        return start_pos, end_pos

    # Try extend left
    new_start_pos = start_pos
    if sequence[start_pos: end_pos].replace("-", "")[:2] == category_heteropolymer:
        while True:
            while new_start_pos - 1 >= 0 and start_pos - new_start_pos < 5 and sequence[new_start_pos - 1] == "-":
                new_start_pos -= 1
            if new_start_pos - 1 >= 0 and sequence[new_start_pos - 1] == category_heteropolymer[-1]:
                new_start_pos -= 1
            else:
                break
            while new_start_pos - 1 >= 0 and start_pos - new_start_pos < 5 and sequence[new_start_pos - 1] == "-":
                new_start_pos -= 1
            if new_start_pos - 1 >= 0 and sequence[new_start_pos - 1] == category_heteropolymer[0]:
                new_start_pos -= 1
            else:
                break
            start_pos = new_start_pos

    # Try extend right
    new_end_pos = end_pos
    if sequence[start_pos: end_pos].replace("-", "")[-2:] == category_heteropolymer:
        while True:
            while new_end_pos + 1 < len(sequence) and new_end_pos - end_pos < 5 and sequence[new_end_pos] == "-":
                new_end_pos += 1
            if new_end_pos + 1 < len(sequence) and sequence[new_end_pos] == category_heteropolymer[0]:
                new_end_pos += 1
            else:
                break
            while new_end_pos + 1 < len(sequence)  and new_end_pos - end_pos < 5 and sequence[new_end_pos] == "-":
                new_end_pos += 1
            if new_end_pos + 1 < len(sequence) and sequence[new_end_pos] == category_heteropolymer[1]:
                new_end_pos += 1
            else:
                break
            end_pos = new_end_pos

    return start_pos, end_pos

def update_abundance_error_dict(g: str, r: str, length: int, category: str, dict_temp_abundance):
    """Update dictionaries storing heteropolymer abundance and sequencing error rates
    """
    if isinstance(length, str):
        length = int(length[:-1])

    # Abundance
    if category not in heteropolymer_abundance[gc_category]:
        category = category[::-1]

    dict_temp_abundance[gc_category][category][0] += length / len(genome_aln.replace("-", "")) * 100
    dict_temp_abundance[gc_category][category][1] += 1


#    heteropolymer_abundance[gc_category][category][-1] += int(length/2) / len(genome_aln.replace("-", "")) * 100
    # Error rate
    nb_errors = 0
    for i, base_genome in enumerate(g):
        base_read = r[i]
        if base_read != base_genome:
            nb_errors += 1
    heteropolymer_error_rate[gc_category][category][0] += nb_errors
    heteropolymer_error_rate[gc_category][category][1] += len(g)
    return dict_temp_abundance


def compute_results():
    """Computes all results, and output as graph and raw file
    """

    # Abundance and error rates for each heteropolymer category

    # 1) Save a plot
    width = 0.2
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
#    L_labels = ['AC-CA', 'AG-GA', 'AT-TA', 'CG-GC', 'CT-TC', 'GT-TG'] # alphabetical sort
    L_labels = ['CT-TC', 'AT-TA', 'GT-TG', 'CG-GC', 'AC-CA', 'AG-GA']

    ## Primary axis: abundance
    abundance_x = np.arange(1, 7, 1)
    offset = -0.25
    list_colors = [COLOR_LOW_GC, COLOR_HIGH_GC, COLOR_HUMAN]
    for gc_cat in ["low GC", "high GC", "human"]:
        color = list_colors[0]
        list_colors = list_colors[1:]
        abundance_y = []
        for het in list_heteropolymer_category:
            abundance_y += [np.mean(heteropolymer_abundance[gc_cat][het])]
        ax1.bar(abundance_x + offset, abundance_y, width=width, color=color, label=gc_cat)
        offset += 0.25
    ## Secondary axis: error rate
    lw_out = 3
    lw_in = 2.5
    marker_size = 8
    error_x = np.arange(1, 7, 1)
    list_colors = [COLOR_LOW_GC, COLOR_HIGH_GC, COLOR_HUMAN]
    for gc_cat in ["low GC", "high GC", "human"]:
        color = list_colors[0]
        list_colors = list_colors[1:]
        errors_y = []
        for het in list_heteropolymer_category:
            if heteropolymer_error_rate[gc_cat][het][1] == 0:
                errors_y += [-1]
            else:
                errors_y += [heteropolymer_error_rate[gc_cat][het][0] / heteropolymer_error_rate[gc_cat][het][1] * 100]
        ax2.plot(error_x, errors_y, color="k", lw=lw_out)
        ax2.plot(error_x, errors_y, "o-", mec="k", ms=marker_size, lw=lw_in, color=color)
    ## Global settings
    plt.xticks(error_x, L_labels)
    ax1.set_ylim(0, 0.05)
    ax1.set_xlabel("Heteropolymer type")
    ax1.set_ylabel("Heteropolymer occurrences in reference genome (%)")
    ax2.set_ylabel("Error rate (%)")

#    ax1.set_ylim(0, 1)
    ax2.set_ylim(2, 17)
    ax2.yaxis.set_minor_locator(MultipleLocator(1))
    ax2.tick_params(which='minor', length=2)

    low_patch = mpatches.Patch(color=COLOR_LOW_GC, label='Low GC bacteria')
    high_patch = mpatches.Patch(color=COLOR_HIGH_GC, label='High GC bacteria')
    human_patch = mpatches.Patch(color=COLOR_HUMAN, label='Human')
    abundance_patch = mpatches.Patch(color="gray", label='Occurrences')
    error_patch = Line2D([0], [0], marker="o", ms = 10, color="gray", lw=4,
                         label='Error rate', mec="k")
    plt.legend(handles=[low_patch, high_patch, human_patch,
                        abundance_patch, error_patch],
               loc=2, ncol=2)
    plt.savefig(OUTPUT_PLOT + "heteropolymer_abundance_error_rate.png")
    plt.close()

    # 2) save raw results in .txt file
    RAW_OUTPUT_FILE = open(OUTPUT_RAW + "heteropolymer_abundance_error_rate.txt", "w")
    for gc_cat in ["low GC", "high GC", "human"]:
        RAW_OUTPUT_FILE.write(f"---\nResults for {gc_cat}\n")
        RAW_OUTPUT_FILE.write("\t".join(["Heteropolymer category", "Abundance", "Error rate"]) + "\n")
        for het in list_heteropolymer_category:
            abundance = np.mean(heteropolymer_abundance[gc_cat][het])
            if heteropolymer_error_rate[gc_cat][het][1] == 0:
                error_rate = -1
            else:
                error_rate = round(heteropolymer_error_rate[gc_cat][het][0] / heteropolymer_error_rate[gc_cat][het][1] * 100, 2)
            het = het + "-" + het[::-1]
            RAW_OUTPUT_FILE.write("\t".join([het, str(abundance), str(error_rate)]) + "\n")
    RAW_OUTPUT_FILE.close()


    # Plots dictionary storing heteropolymer length differences
    #     between genomic expected one, and the actual sequenced one
    list_colors = [COLOR_LOW_GC, COLOR_HIGH_GC, COLOR_HUMAN]
    offset = -0.25
    _, ax = plt.subplots()
    for gc_cat in ["low GC", "high GC", "human"]:
        color = list_colors[0]
        list_colors = list_colors[1:]
        L_to_plot = []
        if diff_heteropol_length_sequenced_dict[gc_cat][4] == {}: # skip if no data recorded yet
            continue
        for heteropolymer_length_genome in diff_heteropol_length_sequenced_dict[gc_cat]:
            L_to_plot += [get_statistics(diff_heteropol_length_sequenced_dict[gc_cat][heteropolymer_length_genome])]
        boxprops = dict(linewidth=4)
        box = ax.bxp(L_to_plot, showfliers=False,
                     positions=np.arange(MIN_HETEROPOLYMER_LENGTH, MAX_HETEROPOLYMER_LENGTH, 2)+offset,
                     patch_artist=True,
                     widths=BOXPLOT_WIDTH,
                     boxprops=boxprops)
        for item in ['boxes', 'fliers', 'whiskers', 'caps']:
            plt.setp(box[item], color=color)
        for patch in box['boxes']:
            patch.set(facecolor=color)
        for item in ['medians']:
            plt.setp(box[item], color="black")
        offset += 0.25
    # Plots expected distribution x = y
    L_x, L_y = [], []
    for i in range(MIN_HETEROPOLYMER_LENGTH, MAX_HETEROPOLYMER_LENGTH, 2):
        L_x += [i - 0.2, i, i+0.2]
        L_y += [i, i, i]
    plt.plot(L_x, L_y, "--", color="black")
    # Legend
    low_patch = mpatches.Patch(color=COLOR_LOW_GC, label='Low GC bacteria')
    high_patch = mpatches.Patch(color=COLOR_HIGH_GC, label='High GC bacteria')
    human_patch = mpatches.Patch(color=COLOR_HUMAN, label='Human')
    plt.legend(handles=[low_patch, high_patch, human_patch],
               title="Species:", loc=9, ncol=3)
    # Details of the plot
    plt.xlim(3.5, 10.5)
    plt.ylim(0.85, 200)
    plt.yscale("log")
    L_ticks = [0.9, 2, 5, 10, 20, 50, 100, 200]
    L_ticks_labels = [0, 2, 5, 10, 20, 50, 100, 200]
    plt.yticks(np.asarray(L_ticks), labels=L_ticks_labels)
    plt.xticks(ticks=np.arange(MIN_HETEROPOLYMER_LENGTH, MAX_HETEROPOLYMER_LENGTH, 2),
               labels=np.arange(MIN_HETEROPOLYMER_LENGTH, MAX_HETEROPOLYMER_LENGTH, 2))
    plt.xlabel("Expected heteropolymer lengths (reference genome)")
    plt.ylabel("Sequenced heteropolymer lengths (reads)")
    plt.savefig(OUTPUT_PLOT + "difference_expected_sequenced_heteropolymer_length.png")
    plt.close()

    # Save raw results in .txt file
    RAW_OUTPUT_FILE = open(OUTPUT_RAW + "difference_expected_sequenced_heteropolymer_length.txt", "w")
    for gc_cat in ["low GC", "high GC", "human"]:
        RAW_OUTPUT_FILE.write(f"---\nResults for {gc_cat}\n")
        RAW_OUTPUT_FILE.write("\t".join(["Heteropolymer length",
                                     "Minimum", "Q1", "Median", "Q3", "Maximum"]) + "\n")
        for heteropolymer_length in range(MIN_HETEROPOLYMER_LENGTH, MAX_HETEROPOLYMER_LENGTH, 2):
            L_to_plot = []
            if diff_heteropol_length_sequenced_dict[gc_cat][4] == {}: # skip if no data recorded yet
                continue
#            for heteropolymer_length_genome in diff_heteropol_length_sequenced_dict[gc_cat]:
            result = get_statistics(diff_heteropol_length_sequenced_dict[gc_cat][heteropolymer_length])
            if result["whislo"] == 0.9: # artefact for plotting log-scaled
                result["whislo"] = 0
            RAW_OUTPUT_FILE.write(f"{heteropolymer_length}\t{result['whislo']}\t{result['q1']}\t{result['med']}\t{result['q3']}\t{result['whishi']}\n")
        RAW_OUTPUT_FILE.write("\n")
    RAW_OUTPUT_FILE.close()


# --------------------------------------------------------------------------------------------------
# Main

if __name__ == "__main__":

    ## Parses arguments
    if len(sys.argv) == 1:
        ALN_EXPL_DIRNAME = "/home/cdelahay/Pipeline/Aln"
        REFERENCE_GENOME_DIRNAME = "/home/cdelahay/Pipeline/Ref"
        OUTPUT_RAW = "/home/cdelahay/Pipeline/Out_raw"
        OUTPUT_PLOT = "/home/cdelahay/Pipeline/Out_graph"
        FILENAME_SPECIES_GC = "/home/cdelahay/Pipeline/species_color_GC.txt"


    else:
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

    MIN_HETEROPOLYMER_LENGTH, MAX_HETEROPOLYMER_LENGTH = 4, 12 # max excluded

    COLOR_LOW_GC = "#4477AA"
    COLOR_HIGH_GC = "#EE6677"
    COLOR_HUMAN = "#CCBB44"

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


    # Width of boxplot
    BOXPLOT_WIDTH = 0.2

    #list_heteropolymer_category = ["AC", "AG", "AT", "CG", "CT", "GT"] # alphabetically sorted
    list_heteropolymer_category = ["CT", "AT", "GT", "CG", "AC", "AG"]
    list_genomic_heteropolymer_lengths = list(np.arange(4, 12, 2)) + ["12+"]

    # Initiate dictionary that will store genomic distribution of heteropolymers
    genomic_distribution_heteropolymer_dict = {}
    for length in list_genomic_heteropolymer_lengths:
        genomic_distribution_heteropolymer_dict[length] = {}
        for category in list_heteropolymer_category:
            genomic_distribution_heteropolymer_dict[length][category] = 0

    # Initiate dictionary that will store, for each genomic heteropolymer's length
    #  the length of sequenced heteropolymer
    diff_heteropol_length_sequenced_dict = initiate_diff_len_dict()


    # Initiate dictionaries that will store genomic (aligned) abundance for each heteropolymer
    #  category
    heteropolymer_abundance = {}
    for gc_category in set(dict_species_gc_category.values()):
        heteropolymer_abundance[gc_category] = {}
        for heteropolymer_category in list_heteropolymer_category:
            heteropolymer_abundance[gc_category][heteropolymer_category] = []

    # Initiate dictionary that will store error rate for each hetoropolymer category
    heteropolymer_error_rate = {}
    for gc_category in set(dict_species_gc_category.values()):
        heteropolymer_error_rate[gc_category] = {}
        for heteropolymer_category in list_heteropolymer_category:
            heteropolymer_error_rate[gc_category][heteropolymer_category] = [0, 0] # nb of errors and total base nb


    # Pattern for regular expression finding of heteropolymers
    pattern_heteropolymer = generate_heteropolymer_pattern()

    species_counter = 0
    for aln_filename in os.listdir(ALN_EXPL_DIRNAME):
        species_counter += 1

        aln_path = ALN_EXPL_DIRNAME + aln_filename
        NB_TOT_ALN = get_total_number_of_lines(aln_path)
        STARTING_TIME = time.time()

        print(aln_filename)
        species_name = aln_filename.split(".txt")[0]
        aln_file = open(aln_path, "r")

        gc_category = dict_species_gc_category[species_name]

        # Get heteropolymer distribution from reference genome
        reference_genome_filename = REFERENCE_GENOME_DIRNAME + species_name + ".fasta"
        reference_genome_length = get_genomic_heteropolymer_distribution()

        # Get alignment from aln_file
        nb_aln_done = 0
        progressing = 0

        dict_temp_heteropol_abundance = {}
        dict_temp_heteropol_abundance[gc_category] = {}
        for het in list_heteropolymer_category:
            dict_temp_heteropol_abundance[gc_category][het] = [0, 0] # sum, nb occurrences


        while True:
            header_aln = aln_file.readline().replace("\n", "")
            if not header_aln:
                break
            genome_aln = aln_file.readline().replace("\n", "")
            read_aln = aln_file.readline().replace("\n", "")

            nb_aln_done += 3
            tmp_progressing = int(nb_aln_done / NB_TOT_ALN * 100)
            if tmp_progressing > progressing and tmp_progressing % 5 == 0:
                progressing = tmp_progressing
                display_progressing_bar(progressing, time.time())
                compute_results()

            if nb_aln_done > NB_TOT_ALN:
                break

            # Compute sequencing errors associated with genomic heteropolymers
            start_all, end_all = -1, -1
            for res in re.finditer(pattern_heteropolymer, genome_aln):
                start, end = res.span()
                if start in range(start_all, end_all): # this heteropolymer already has been analysed
                    continue
                category_heteropolymer = res.group()[:2]
                start_genome, end_genome = try_extend_genomic_heteropolymer(start, end, genome_aln)
                heteropolymer_genome = genome_aln[start_genome: end_genome]
                heteropolymer_genome_length = len(heteropolymer_genome.replace("-", ""))
                if heteropolymer_genome_length >= 12:
                    heteropolymer_genome_length = "12+"
                start_read, end_read = try_extend_read_heteropolymer(start, end, read_aln)
                start_all = min(start_genome, start_read)
                end_all = max(end_genome, end_read)

                # Update dictionary storing heteropolymer length differences
                #     between genomic expected one, and the actual sequenced one
                update_diff_len_dict(genome_aln[start_genome: end_genome],
                                     read_aln[start_read: end_read],
                                     heteropolymer_genome_length,
                                     category_heteropolymer)

                # Update dictionaries storing heteropolymer abundance and error rates
                dict_temp_heteropol_abundance = update_abundance_error_dict(genome_aln[start_all: end_all],
                                                                           read_aln[start_all: end_all],
                                                                           heteropolymer_genome_length,
                                                                           category_heteropolymer,
                                                                           dict_temp_heteropol_abundance)

        for het in dict_temp_heteropol_abundance[gc_category]:
            sum_het_abundance, occ_het_abundance = dict_temp_heteropol_abundance[gc_category][het]
            heteropolymer_abundance[gc_category][het] += [sum_het_abundance / occ_het_abundance]
        aln_file.close()
        sys.stdout.write("\n")


    compute_results()
