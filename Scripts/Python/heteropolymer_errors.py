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

try:
    import seaborn as sns; sns.set_theme()
except ImportError:
    print("Package seaborn is required, please install it")
    sys.exit(1)


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


def generate_heteropolymer_pattern():
    """Generates a regular expression to find heteropolymers
    i.e. repetition of groups of 2 nucleotides
    """
    L_bases = ["A", "C", "G", "T"]
    pattern = ""
    for b1 in L_bases:
        for b2 in L_bases:
            if b1 == b2:
                continue
            pattern += f"(?P<{b1+b2}>({b1+b2}){{2,}})|"
    pattern = pattern[:-1]
    return pattern


def initiate_diff_len_dict() -> dict:
    """Add a 'sub-dictionary' for the current species, which will, for each species group and
    for each genomic heteropolymer length, store the sequenced heteropolymer length (i.e. in the read)
    """
    dictionary = {}
    for group_name in L_groups:
        dictionary[group_name] = {}
        for error_length in L_heteropolymer_lengths:
            dictionary[group_name][error_length] = {}
    return dictionary


def get_total_number_of_lines(filename: str) -> int:
    """Returns number of lines of filename given as input
    """
    with open(filename, "r") as file:
        return sum(1 for _ in file)

#
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


def get_genomic_heteropolymer_distribution(group_name):
    """Compute heteropolymer distribution in the current reference genome, and
    updates the associated dictionary
    Also returns length of genome
    """
    L_labels = ['AC-CA', 'AG-GA', 'AT-TA', 'CG-GC', 'CT-TC', 'GT-TG']

    # Colors depending on heteropolymer length
    list_colors = sns.diverging_palette(240, 10, n=len(L_heteropolymer_lengths), as_cmap=False)
    list_colors_hex = list_colors.as_hex() # get hexadecimal code for each color
    dict_colors = {}
    for i, length in enumerate(L_heteropolymer_lengths):
        dict_colors[length] = list_colors_hex[i]

    # Add "space" to store results
    for length in dict_genomic_distrib[group_name]:
        for heteropolymer in L_heteropolymers:
            dict_genomic_distrib[group_name][length][heteropolymer] = 0


    # Get reference genome and its length
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

    # Get genomic distribution, with regex applied to genome
    for res in re.finditer(pattern, genome):
        start, end = res.span()
        heteropolymer_length = end - start
        if heteropolymer_length >= MAX_HETEROPOLYMER_LENGTH:
            heteropolymer_length = f"{MAX_HETEROPOLYMER_LENGTH}+"
        base = res.group()[:2]
        # use canonical form:
        if base[::-1] < base:
            base = base[::-1]
        dict_genomic_distrib[group_name][heteropolymer_length][base] += 1

    # Prepare data for plot
    for group_name in L_groups:
        dict_to_plot = {}
        if dict_genomic_distrib[group_name][MIN_HETEROPOLYMER_LENGTH] == {}:
            continue
        for length in dict_genomic_distrib[group_name]:
            dict_to_plot[length] = []
            for heteropolymer in L_heteropolymers:
                occurrence = round(dict_genomic_distrib[group_name][length][heteropolymer], 2)
                dict_to_plot[length] += [occurrence]


        # Plot results:
        fig, ax = plt.subplots(constrained_layout=True)
        r = np.arange(len(L_heteropolymers))
        for length in L_heteropolymer_lengths:
            r = [x + BAR_WIDTH for x in r]
            ax.bar(r, dict_to_plot[length],
                   color=dict_colors[length],
                   width=BAR_WIDTH,
                   edgecolor="white", label=length)
        ax.set_yscale('log')
        ax.set_xlabel("Heteropolymer type", fontweight="bold")
        ax.set_ylabel("Occurrences", fontweight="bold")
        plt.xticks([r + 3*BAR_WIDTH for r in range(len(L_labels))], L_labels)
        plt.ylim(1e2, 5*1e8)
        ax.legend(ncol=5, title="Heteropolymer length", loc=1)
        plt.savefig(OUTPUT_PLOT + f"heteropolymer_genomic_distribution_{group_name}.png")
        plt.close()

    # Save raw results in .txt file
    RAW_OUTPUT_FILE = open(OUTPUT_RAW + "heteropolymer_genomic_distribution.txt", "w")
    for group_name in L_groups:
        RAW_OUTPUT_FILE.write(f"{group_name}\n")
        if dict_genomic_distrib[group_name][MIN_HETEROPOLYMER_LENGTH] == {}:
            continue
        RAW_OUTPUT_FILE.write("\t".join(["Heteropolymer length"] + L_labels) + "\n")
        for heteropolymer_length in dict_genomic_distrib[group_name]:
            RAW_OUTPUT_FILE.write(str(heteropolymer_length))
            for heteropolymer in L_heteropolymers:
                occurrence = np.median(dict_genomic_distrib[group_name][heteropolymer_length][heteropolymer])
                RAW_OUTPUT_FILE.write("\t" + str(occurrence))
            RAW_OUTPUT_FILE.write("\n")
        RAW_OUTPUT_FILE.write("\n")
    RAW_OUTPUT_FILE.close()

    return reference_genome_length


def update_diff_len_dict(g: str, r: str, heteropol_len_genome: int, heteropol_category: str, group_name: str):
    """Updates dictionary that stores heteropolymer length differences (genomic and sequenced)
    """
    g = g.replace("-", "")
    r = r.replace("-", "")
    if heteropol_len_genome not in dict_diff_length_sequenced[group_name]:
        return
    if heteropol_category not in r:
        heteropol_len_read = 0
    else:
        potential_heteropol_len_read = []
        pattern_tmp = "(" + heteropol_category + "){0,}"
        for match in re.finditer(pattern_tmp, r):
            match_length = match.span()[1] - match.span()[0]
            potential_heteropol_len_read += [match_length]
        heteropol_len_read = max(potential_heteropol_len_read)
    if heteropol_len_read not in dict_diff_length_sequenced[group_name][heteropol_len_genome]:
        dict_diff_length_sequenced[group_name][heteropol_len_genome][heteropol_len_read] = 0
    dict_diff_length_sequenced[group_name][heteropol_len_genome][heteropol_len_read] += 1


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


def update_abundance_error_dict(g: str, r: str, length: int, het_category: str):
    """Update dictionaries storing heteropolymer abundance and sequencing error rates
    """
    if isinstance(length, str):
        length = int(length[:-1])

    # Abundance
    if het_category not in dict_abundance[group_name]:
        het_category = het_category[::-1]

    dict_tmp_abundance[group_name][het_category][0] += length / len(genome_aln.replace("-", "")) * 100
    dict_tmp_abundance[group_name][het_category][1] += 1

    # Error rate
    nb_errors = 0
    for i, base_genome in enumerate(g):
        base_read = r[i]
        if base_read != base_genome:
            nb_errors += 1
    dict_error_rate[group_name][het_category][0] += nb_errors
    dict_error_rate[group_name][het_category][1] += len(g)


def compute_results():
    """Computes all results, and output as graph and raw file
    """

    # --- Abundance and error rates for each heteropolymer category ---

    # 1/ Save as plot
    width = 0.2
    L_labels = ['AC-CA', 'AG-GA', 'AT-TA', 'CG-GC', 'CT-TC', 'GT-TG'] # alphabetical sort
#    L_labels = ['CT-TC', 'AT-TA', 'GT-TG', 'CG-GC', 'AC-CA', 'AG-GA']

    # Primary axis = abundance ; secondary axis = error rate
    abundance_x = np.arange(1, len(L_heteropolymers)+1, 1)
    error_x = np.arange(1, len(L_heteropolymers)+1, 1)

    for group_name in L_groups:
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        L_patches = []
        abundance_y = []
        errors_y = []
        for heteropolymer in L_heteropolymers:
            abundance_y += [np.mean(dict_abundance[group_name][heteropolymer])]
            if dict_error_rate[group_name][heteropolymer][1] == 0:
                errors_y += [-1]
            else:
                errors_y += [dict_error_rate[group_name][heteropolymer][0] / \
                             dict_error_rate[group_name][heteropolymer][1] * 100]
        ax2.plot(error_x, errors_y, color="k")
        ax2.plot(error_x, errors_y, "o-", mec="k")
        ax1.bar(abundance_x, abundance_y, width=width, label=group_name)
        patch = mpatches.Patch(label=group_name)
        abundance_patch = mpatches.Patch(color="gray", label='Occurrences')
        error_patch = Line2D([0], [0], marker="o", color="gray", label='Error rate', mec="k")
        L_patches = [patch, abundance_patch, error_patch]

        ## Global settings
        plt.xticks(error_x, L_labels)
        ax1.set_xlabel("Heteropolymer type")
        ax1.set_ylabel("Heteropolymer occurrences in reference genome (%)")
        ax2.set_ylabel("Error rate (%)")

        ax2.yaxis.set_minor_locator(MultipleLocator(1))
        ax2.tick_params(which='minor', length=2)

        plt.legend(handles=L_patches)
        plt.savefig(OUTPUT_PLOT + f"heteropolymer_abundance_error_rate_{group_name}.png")
        plt.close()


    # 2/ save raw results in .txt file
    RAW_OUTPUT_FILE = open(OUTPUT_RAW + "heteropolymer_abundance_error_rate.txt", "w")
    for group_name in L_groups:
        RAW_OUTPUT_FILE.write(f"\n{group_name}\n")
        RAW_OUTPUT_FILE.write("\t".join(["Heteropolymer category", "Abundance", "Error rate"]) + "\n")
        for heteropolymer in L_heteropolymers:
            abundance = np.mean(dict_abundance[group_name][heteropolymer])
            if dict_error_rate[group_name][heteropolymer][1] == 0:
                error_rate = -1
            else:
                error_rate = round(dict_error_rate[group_name][heteropolymer][0] / \
                                   dict_error_rate[group_name][heteropolymer][1] * 100, 2)
            heteropolymer = heteropolymer + "-" + heteropolymer[::-1]
            RAW_OUTPUT_FILE.write("\t".join([heteropolymer, str(abundance), str(error_rate)]) + "\n")
    RAW_OUTPUT_FILE.close()


    # --- Heteropolymer length differences ---
    #     between genomic expected one, and the actual sequenced one

    for group_name in L_groups:
        _, ax = plt.subplots()
        L_patches = []
        color = dict_species_group_color[group_name]
        L_to_plot = []
        if dict_diff_length_sequenced[group_name][MIN_HETEROPOLYMER_LENGTH] == {}: # skip if no data recorded yet
            continue
        for heteropolymer_length_genome in dict_diff_length_sequenced[group_name]:
            if dict_diff_length_sequenced[group_name][heteropolymer_length_genome] == {}:
                break
            L_to_plot += [get_statistics(dict_diff_length_sequenced[group_name][heteropolymer_length_genome])]
        boxprops = dict(linewidth=4)
        box = ax.bxp(L_to_plot, showfliers=False,
                     positions=np.arange(MIN_HETEROPOLYMER_LENGTH, MAX_HETEROPOLYMER_LENGTH+1, 2)[:len(L_to_plot)],
                     patch_artist=True, widths=BOXPLOT_WIDTH, boxprops=boxprops)
        for item in ['boxes', 'fliers', 'whiskers', 'caps']:
            plt.setp(box[item], color=color)
        for patch in box['boxes']:
            patch.set(facecolor=color)
        for item in ['medians']:
            plt.setp(box[item], color="black")
        # Plots expected distribution x = y
        L_x, L_y = [], []
        for i in range(MIN_HETEROPOLYMER_LENGTH, MAX_HETEROPOLYMER_LENGTH+1, 2):
            L_x += [i - 0.2, i, i+0.2]
            L_y += [i, i, i]
        plt.plot(L_x, L_y, "--", color="black")
        # Details of the plot
        plt.yscale("log")
        plt.xticks(ticks=np.arange(MIN_HETEROPOLYMER_LENGTH, MAX_HETEROPOLYMER_LENGTH+1, 2),
                   labels=L_heteropolymer_lengths)
        plt.xlabel("Expected heteropolymer lengths (reference genome)")
        plt.ylabel("Sequenced heteropolymer lengths (reads)")
        plt.savefig(OUTPUT_PLOT + f"difference_expected_sequenced_heteropolymer_length_{group_name}.png")
        plt.close()


    # Save raw results in .txt file
    RAW_OUTPUT_FILE = open(OUTPUT_RAW + "difference_expected_sequenced_heteropolymer_length.txt", "w")
    for group_name in L_groups:
        RAW_OUTPUT_FILE.write(f"\n{group_name}\n")
        RAW_OUTPUT_FILE.write("\t".join(["Heteropolymer length",
                                     "Minimum", "Q1", "Median", "Q3", "Maximum"]) + "\n")
        for heteropolymer_length in L_heteropolymer_lengths:
            L_to_plot = []
            if dict_diff_length_sequenced[group_name][heteropolymer_length] == {}: # skip if no data recorded
                continue
            result = get_statistics(dict_diff_length_sequenced[group_name][heteropolymer_length])
            if result["whislo"] == 0.9: # artefact for plotting log-scaled
                result["whislo"] = 0
            RAW_OUTPUT_FILE.write(f"{heteropolymer_length}\t{result['whislo']}\t{result['q1']}\t{result['med']}\t{result['q3']}\t{result['whishi']}\n")
    RAW_OUTPUT_FILE.close()


# --------------------------------------------------------------------------------------------------
# Main

if __name__ == "__main__":

    ## Parses arguments
    NUMBER_EXPECTED_ARGUMENTS = 8
    if len(sys.argv) != NUMBER_EXPECTED_ARGUMENTS + 1:
        print(f"ERROR: Wrong number of arguments: {NUMBER_EXPECTED_ARGUMENTS} expected but {len(sys.argv)-1} given.")
        sys.exit(2)
    ALN_EXPL_DIRNAME, REFERENCE_GENOME_DIRNAME, OUTPUT_RAW, OUTPUT_PLOT, FILE_SPECIES_GC, FILE_SPECIES_GROUP, MIN_HETEROPOLYMER_LENGTH, MAX_HETEROPOLYMER_LENGTH = sys.argv[1:]

    MIN_HETEROPOLYMER_LENGTH = int(MIN_HETEROPOLYMER_LENGTH)
    MAX_HETEROPOLYMER_LENGTH = int(MAX_HETEROPOLYMER_LENGTH)

    # --- Parameters ---

    # Length of heteropolymers
    L_heteropolymer_lengths = list(np.arange(MIN_HETEROPOLYMER_LENGTH, MAX_HETEROPOLYMER_LENGTH, 2)) + [f"{MAX_HETEROPOLYMER_LENGTH}+"]

    # Get groups of species
    dict_species_group, dict_species_group_color = get_groups(FILE_SPECIES_GROUP)
    L_groups = sorted(set(dict_species_group.values()))

    # Width of boxplot
    BOXPLOT_WIDTH = 0.2
    BAR_WIDTH = 0.1

    L_heteropolymers = ["AC", "AG", "AT", "CG", "CT", "GT"]

    # Pattern for regular expression finding of heteropolymers
    pattern = generate_heteropolymer_pattern()

    # --- Initiate dictionaries that will store results ---

    # Store genomic distribution of heteropolymers, for each heteropolymer
    dict_genomic_distrib = {}
    for group_name in L_groups:
        dict_genomic_distrib[group_name] = {}
        for length in L_heteropolymer_lengths:
            dict_genomic_distrib[group_name][length] = {}

    # Store, for each genomic heteropolymer's length, the length of sequenced heteropolymer
    dict_diff_length_sequenced = initiate_diff_len_dict()

    # Store genomic (aligned) abundance for each heteropolymer category
    dict_abundance = {}
    for group_name in L_groups:
        dict_abundance[group_name] = {}
        for heteropolymer_category in L_heteropolymers:
            dict_abundance[group_name][heteropolymer_category] = []

    # Store error rate for each hetoropolymer category
    #    dictionary[group][heteropolymer] = # [number of errors, total base count]
    dict_error_rate = {}
    for group_name in L_groups:
        dict_error_rate[group_name] = {}
        for heteropolymer_category in L_heteropolymers:
            dict_error_rate[group_name][heteropolymer_category] = [0, 0]


    # --- Analysis ---
#
    for aln_filename in os.listdir(ALN_EXPL_DIRNAME):

        aln_path = ALN_EXPL_DIRNAME + aln_filename

        # Get total number of alignments in this file
        #     divided by 3 as alignments are represented as triplets of lines
        NB_TOT_ALN = int(get_total_number_of_lines(aln_path) / 3)
        STARTING_TIME = time.time()

        species_name = aln_filename.split(".txt")[0]
        print("  ", species_name)
        aln_file = open(aln_path, "r")

        group_name = dict_species_group[species_name]

        # Get heteropolymer distribution from reference genome
        reference_genome_filename = REFERENCE_GENOME_DIRNAME + species_name + ".fasta"
        reference_genome_length = get_genomic_heteropolymer_distribution(group_name)

        # Get alignment from aln_file
        nb_aln_done = 0
        progressing = 0

        # Use temporary dictionary to store abundance
        #   dict_tmp[group][heterpolymer] = [sum of occurrences, count of occurrences]
        dict_tmp_abundance = {}
        dict_tmp_abundance[group_name] = {}
        for heteropolymer in L_heteropolymers:
            dict_tmp_abundance[group_name][heteropolymer] = [0, 0]

        # Browse alignments
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


            if nb_aln_done > 100: ## TODO remove
                break


            # Compute sequencing errors associated with genomic heteropolymers
            start_all, end_all = -1, -1
            for res in re.finditer(pattern, genome_aln):
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
                                     category_heteropolymer,
                                     group_name)

                # Update dictionaries storing heteropolymer abundance and error rates
                update_abundance_error_dict(genome_aln[start_all: end_all],
                                            read_aln[start_all: end_all],
                                            heteropolymer_genome_length,
                                            category_heteropolymer)

        for heteropolymer in dict_tmp_abundance[group_name]:
            sum_abundance, occurence_abundance = dict_tmp_abundance[group_name][heteropolymer]
            dict_abundance[group_name][heteropolymer] += [sum_abundance / occurence_abundance]
        aln_file.close()
        sys.stdout.write("\n")


    compute_results()
