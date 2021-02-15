#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Analyses trinucleotides, i.e. repetitions of trio of bases.
Output (plot and raw), depending on number of repetitions and trinucleotide category:
    - distribution in the genome
    - sequencing accuracy
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

import matplotlib.ticker as ticker




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



def get_regex_pattern():
    """Returns a regex pattern to find trinucleotides
    """
    pattern = ""
    for base1 in LIST_BASES:
        for base2 in LIST_BASES:
            for base3 in LIST_BASES:
                #trinucleotide = base1 + base2 + base3 # not taking into account potential gaps in aligned genome…
                if base1 == base2 == base3: # Omit homopolymers
                    continue
                trinucleotide = "-*".join([base1, base2, base3, ""]) # allow gaps within pattern in aligned genome
                pattern += f"({trinucleotide}){{2,}}|"
    pattern = pattern[:-1]
    return pattern

def convert_trinucleotide_to_SW_pattern(trinucleotide: str):
    """Converts a trinucleotide to a pattern made of S and W
    S (strong) for bases C or G
    W (weak) for bases A or T
    """
    sw_pattern = ""
    for base in trinucleotide:
        if base in ["C", "G"]:
            sw_pattern += "S"
        else:
            sw_pattern += "W"
    return sw_pattern

def convert_sw_pattern_to_global_pattern(sw_pattern: str):
    """Converts a SW pattern (three letters in [W, S] to  "S-only", "Mostly-S", "W-only", and
    "Mosly-W" global pattern
    """
    count_W = sw_pattern.count("W")
    if count_W == 0:
        sw_global_pattern = "S only"
    elif count_W == 1:
        sw_global_pattern = "mainly S"
    elif count_W == 2:
        sw_global_pattern = "mainly W"
    else:
        sw_global_pattern = "W only"
    return sw_global_pattern

def initiate_result_dict():
    """Initialize dictionary that will store every information needed to plot both
    distribution of trinucleotides depending on category and length
    and sequencing accuracy
    For each gc category, each trinucleotide type and each length, store two lists:
        - error rate for each species
        - occurrences for each species
    """
    dictionary = {}
    for base1 in LIST_BASES:
        for base2 in LIST_BASES:
            for base3 in LIST_BASES:
                trinucleotide = base1 + base2 + base3
                sw_pattern = convert_trinucleotide_to_SW_pattern(trinucleotide)
                sw_pattern = convert_sw_pattern_to_global_pattern(sw_pattern)
                for gc_category in ["low GC", "high GC", "human"]:
                    if gc_category not in dictionary:
                        dictionary[gc_category] = {}
                    if sw_pattern not in dictionary[gc_category]:
                        dictionary[gc_category][sw_pattern] = {}
                        dictionary[gc_category][sw_pattern]["4+"] = [[], []]
                        for l in range(2, 4):
                            dictionary[gc_category][sw_pattern][l] = [[], []]
    return dictionary


def initiate_temp_dict():
    """Initialize TEMPORARY dictionary that will store every information needed to plot both
    distribution of trinucleotides depending on category and length
    and sequencing accuracy
    For each trinucleotide type and each length, store a list of 3 integers:
        - number of erroneously sequenced bases
        - number of sequenced bases
        - occurrences of such trinucleotide
    """
    dictionary = {}
    for base1 in LIST_BASES:
        for base2 in LIST_BASES:
            for base3 in LIST_BASES:
                trinucleotide = base1 + base2 + base3
                sw_pattern = convert_trinucleotide_to_SW_pattern(trinucleotide)
                sw_pattern = convert_sw_pattern_to_global_pattern(sw_pattern)
                dictionary[sw_pattern] = {}
                dictionary[sw_pattern]["4+"] = [0, 0, 0]
                for l in range(2, 4):
                    dictionary[sw_pattern][l] = [0, 0, 0]
    return dictionary

def is_homopolymer(word: str) -> bool:
    """Returns True if the given word (trinucleotide) is an homopolymer, i.e. a repetition of
    same base. Else returns False.
    Homopolymers are excluded here, as they have a dedicated study.
    """
    if word[0] * len(word) == word:
        return True
    return False


def try_extend_read_trinucleotide(start_pos, end_pos, sequence):
    """ If sequence[start: end] start/end by '-' or trinucleotide category, extend the
    start and end positions (delimiting sequenced trinucleotide)
    + apply constraints on maximum gap length
    """
    if sequence[start_pos: end_pos].replace("-", "") == "":
        return start_pos, end_pos

    # Try extend left
    new_start_pos = start_pos
    if sequence[start_pos: end_pos].replace("-", "")[:3] == category:
        while True:
            while new_start_pos - 1 >= 0 and start_pos - new_start_pos < 5 and sequence[new_start_pos - 1] == "-":
                new_start_pos -= 1
            if new_start_pos - 1 >= 0 and sequence[new_start_pos - 1] == category[-1]:
                new_start_pos -= 1
            else:
                break
            while new_start_pos - 1 >= 0 and start_pos - new_start_pos < 5 and sequence[new_start_pos - 1] == "-":
                new_start_pos -= 1
            if new_start_pos - 1 >= 0 and sequence[new_start_pos - 1] == category[-2]:
                new_start_pos -= 1
            else:
                break
            while new_start_pos - 1 >= 0 and start_pos - new_start_pos < 5 and sequence[new_start_pos - 1] == "-":
                new_start_pos -= 1
            if new_start_pos - 1 >= 0 and sequence[new_start_pos - 1] == category[0]:
                new_start_pos -= 1
            else:
                break
            start_pos = new_start_pos

    # Try extend right
    new_end_pos = end_pos
    if sequence[start_pos: end_pos].replace("-", "")[-3:] == category:
        while True:
            while new_end_pos + 1 < len(sequence) and new_end_pos - end_pos < 5 and sequence[new_end_pos] == "-":
                new_end_pos += 1
            if new_end_pos + 1 < len(sequence) and sequence[new_end_pos] == category[0]:
                new_end_pos += 1
            else:
                break
            while new_end_pos + 1 < len(sequence) and new_end_pos - end_pos < 5 and sequence[new_end_pos] == "-":
                new_end_pos += 1
            if new_end_pos + 1 < len(sequence) and sequence[new_end_pos] == category[1]:
                new_end_pos += 1
            else:
                break
            while new_end_pos + 1 < len(sequence)  and new_end_pos - end_pos < 5 and sequence[new_end_pos] == "-":
                new_end_pos += 1
            if new_end_pos + 1 < len(sequence) and sequence[new_end_pos] == category[2]:
                new_end_pos += 1
            else:
                break
            end_pos = new_end_pos

    return start_pos, end_pos


def update_results(g: str, r: str):
    """Updates dictionary storing abundance and error rate for each trinucleotide category
    """
    sw_pattern = convert_trinucleotide_to_SW_pattern(category)
    sw_pattern = convert_sw_pattern_to_global_pattern(sw_pattern)
    # Store abundance
    trinucleotide_occ_acc_temp_dict[sw_pattern][trinucleotide_genome_length][2] += 1
    # Store number of erroneous bases
    nb_errors = 0
    for i, base_genome in enumerate(g):
        base_read = r[i]
        if base_read != base_genome:
            nb_errors += 1
    trinucleotide_occ_acc_temp_dict[sw_pattern][trinucleotide_genome_length][0] += nb_errors
    trinucleotide_occ_acc_temp_dict[sw_pattern][trinucleotide_genome_length][1] += len(g)


def compute_results():
    """Compute graph and raw output for current analysis results
    """

    dict_color_trinucl = {"S only": '#8bbeda',
                          "mainly S": '#b2df8a',
                          "mainly W": '#fdbf6f',
                          "W only": '#ea484b'}
    patch_S_only = mpatches.Patch(color=dict_color_trinucl["S only"], label='S only')
    patch_main_S = mpatches.Patch(color=dict_color_trinucl["mainly S"], label='mainly S')
    patch_main_W = mpatches.Patch(color=dict_color_trinucl["mainly W"], label='mainly W')
    patch_W_only = mpatches.Patch(color=dict_color_trinucl["W only"], label='W only')

    width = 0.7
    list_offset = [-0.5, 0.5, 0]

    # Abundance of trinucleotides
    fig, ax = plt.subplots()
    fig.set_dpi(300)
    for i, gc_category in enumerate(trinucleotide_occurrences_accuracy_dict):
        offset = list_offset[i]
        start_pos = np.arange(10, 40, 10) + offset
        for trinucleotide in ["S only", "mainly S", "mainly W", "W only"]:
            list_abundance = []
            for length in [2, 3, "4+"]:
                list_abundance += [trinucleotide_occurrences_accuracy_dict[gc_category][trinucleotide][length][1]]
            if gc_category == "human":
                list_abundance = [np.mean(elt) for elt in list_abundance]
                plt.plot(start_pos, list_abundance, "x", color="k", markersize=10)
            else:
                bplot = ax.boxplot(list_abundance, widths=width,
                                   positions=start_pos,
                                   patch_artist=True,
                                   zorder=-1)
                colors = [dict_color_trinucl[trinucleotide]]*5

                for patch, color in zip(bplot['boxes'], colors):
                    patch.set_facecolor(color)
                    patch.set_edgecolor("black")
                    patch.set_color(color)
                [item.set_color('k') for item in bplot['medians']] # set median line black

            start_pos = start_pos + 3 * width

    ax.legend(handles=[patch_S_only, patch_main_S, patch_main_W, patch_W_only],
              title="Trinucleotides types:", ncol=8)

    plt.yscale("log")
    plt.xlabel("Number of repetitions")
    plt.ylabel("Occurrences ratio of sequenced genomic trinucleotides (%)")
    xtick_positions = [12, 22, 32]
    ax.set_xticks(xtick_positions)
    ax.set_xticklabels([2, 3, "4+"])

    # Set y ticks label as "normal" notation, not scientific one
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y)))

    plt.xlim(8, 38)
    plt.savefig(OUTPUT_PLOT + "trinucleotide_occurrences.png")
    plt.close()


    # Sequencing accuracy of trinucleotides
    fig, ax = plt.subplots()
    fig.set_dpi(300)
    for i, gc_category in enumerate(trinucleotide_occurrences_accuracy_dict):
        offset = list_offset[i]
        start_pos = np.arange(10, 40, 10) + offset
        for trinucleotide in ["S only", "mainly S", "mainly W", "W only"]:
            list_accuracy = []
            for length in [2, 3, "4+"]:
                list_accuracy += [trinucleotide_occurrences_accuracy_dict[gc_category][trinucleotide][length][0]]
            if gc_category == "human":
                list_accuracy = [np.mean(elt) for elt in list_accuracy]
                plt.plot(start_pos, list_accuracy, "x", color="k", markersize=10)
            else:
                bplot = ax.boxplot(list_accuracy, widths=width,
                               positions=start_pos,
                               patch_artist=True,
                               zorder=-1)
                colors = [dict_color_trinucl[trinucleotide]]*5
                for patch, color in zip(bplot['boxes'], colors):
                    patch.set_facecolor(color)
                    patch.set_edgecolor("black")
                    patch.set_color(color)
                [item.set_color('k') for item in bplot['medians']] # set median line black

            start_pos = start_pos + 3 * width

    ax.legend(handles=[patch_S_only, patch_main_S, patch_main_W, patch_W_only],
              title="Trinucleotides types:", ncol=8)

    plt.xlabel("Number of repetitions")
    plt.ylabel("Correctly sequenced trinucleotides (%)")
    xtick_positions = [12, 22, 32]
    ax.set_xticks(xtick_positions)
    ax.set_xticklabels([2, 3, "4+"])
    plt.xlim(8, 38)
    plt.savefig(OUTPUT_PLOT + "trinucleotide_accuracy.png")
    plt.close()


    # Saving results in file
    dict_trinucl_occ = {}
    dict_trinucl_accuracy = {}
    for gc_category in trinucleotide_occurrences_accuracy_dict:
        dict_trinucl_occ[gc_category] = {}
        dict_trinucl_accuracy[gc_category] = {}
        for trinucleotide_category in trinucleotide_occurrences_accuracy_dict[gc_category]:
            dict_trinucl_occ[gc_category][trinucleotide_category] = {}
            dict_trinucl_accuracy[gc_category][trinucleotide_category] = {}
            for length in trinucleotide_occurrences_accuracy_dict[gc_category][trinucleotide_category]:
                list_accuracy, list_occ = trinucleotide_occurrences_accuracy_dict[gc_category][trinucleotide_category][length]
                dict_trinucl_occ[gc_category][trinucleotide_category][length] = list_occ
                dict_trinucl_accuracy[gc_category][trinucleotide_category][length] = list_accuracy


    output = open(OUTPUT_RAW + "trinucleotides_occurrence_accuracy.txt", "w")
    output.write("Occurrences\n")
    output.write(f"{dict_trinucl_occ}\n\n")
    output.write("Accuracy\n")
    output.write(f"{dict_trinucl_accuracy}\n")
    output.close()



# --------------------------------------------------------------------------------------------------
# Main

if __name__ == "__main__":


    ## Parses arguments
    if len(sys.argv) != 5:
        print(f"ERROR: Wrong number of arguments: 4 expected but {len(sys.argv)-1} given.")
        sys.exit(2)
    ALN_EXPL_DIRNAME = sys.argv[1]
    OUTPUT_RAW = sys.argv[2]
    OUTPUT_PLOT = sys.argv[3]
    FILENAME_SPECIES_GC = sys.argv[4]

    if ALN_EXPL_DIRNAME[-1] != "/":
        ALN_EXPL_DIRNAME += "/"
    if OUTPUT_PLOT[-1] != "/":
        OUTPUT_PLOT += "/"
    if OUTPUT_RAW[-1] != "/":
        OUTPUT_RAW += "/"


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

    # Some global variables
    LIST_BASES = ["A", "C", "G", "T"]
    PATTERN = get_regex_pattern() # regex pattern for searching for trinucleotides
    COLOR_LOW_GC = "#4477AA"
    COLOR_HIGH_GC = "#EE6677"
    COLOR_HUMAN = "#CCBB44"


    trinucleotide_occurrences_accuracy_dict = initiate_result_dict() # store results

    # Analyse alignment file
    species_counter = 0
    for aln_filename in os.listdir(ALN_EXPL_DIRNAME):
        species_counter += 1

        aln_path = ALN_EXPL_DIRNAME + aln_filename
        NB_TOT_ALN = get_total_number_of_lines(aln_path) / 3
        STARTING_TIME = time.time()

        print(aln_filename)

        species_name = aln_filename.split(".txt")[0]
        gc_category = dict_species_gc_category[species_name]
        aln_file = open(aln_path, "r")

        trinucleotide_occ_acc_temp_dict = initiate_temp_dict() # Temporary result for each species

        # Get alignment from aln_file
        nb_aln_done = 0
        progressing = 0
        total_nb_sequenced_bases = 0 # To normalise occurrences of trinucleotids
        while True:
            header_aln = aln_file.readline().replace("\n", "")
            if not header_aln:
                break
            genome_aln = aln_file.readline().replace("\n", "")
            read_aln = aln_file.readline().replace("\n", "")

            total_nb_sequenced_bases += len(genome_aln.replace("-", ""))

            nb_aln_done += 1
            tmp_progressing = int(nb_aln_done / NB_TOT_ALN * 100)
            if tmp_progressing > progressing and tmp_progressing % 1 == 0:
                progressing = tmp_progressing
                display_progressing_bar(progressing, time.time())

            #if nb_aln_done >= 100000: # bound number of alignments analysed
            #    break


            # Compute sequencing errors associated with genomic trinucleotides
            start_all, end_all = -1, -1
            for res in re.finditer(PATTERN, genome_aln):
                start, end = res.span()
                if start in range(start_all, end_all): # this trinucleotide has already been analysed
                    continue
                while genome_aln[end - 1] == "-":
                    end -= 1
                category = res.group().replace("-", "")[:3]

                if is_homopolymer(category):
                    continue

                start_genome, end_genome = start, end
                trinucleotide_genome = genome_aln[start_genome: end_genome]


                # length here is the number of repetitions
                trinucleotide_genome_length = len(trinucleotide_genome.replace("-", "")) / 3
                if trinucleotide_genome_length >= 4:
                    trinucleotide_genome_length = "4+"
                start_read, end_read = try_extend_read_trinucleotide(start, end, read_aln)
                start_all = min(start_genome, start_read)
                end_all = max(end_genome, end_read)

                if start != start_genome or end != end_genome:
                    print(genome_aln[start: end])
                    print(category)
                    print(trinucleotide_genome)
                    input("?")


                # Update dictionary storing results
                update_results(genome_aln[start_all: end_all],
                               read_aln[start_all: end_all])


        # Merge temporary result to final dictionary storing results
        threshold = 100
        for sw_pattern in trinucleotide_occ_acc_temp_dict:
            for length in trinucleotide_occ_acc_temp_dict[sw_pattern]:
                nb_err, nb_tot, nb_occ = trinucleotide_occ_acc_temp_dict[sw_pattern][length]
                if nb_occ > threshold: #!= 0:
                    ratio_occ = nb_occ / total_nb_sequenced_bases * 100
                    accuracy_rate = round(100 - (nb_err / nb_tot * 100), 2)
                    trinucleotide_occurrences_accuracy_dict[gc_category][sw_pattern][length][0] += [accuracy_rate]
                    trinucleotide_occurrences_accuracy_dict[gc_category][sw_pattern][length][1] += [ratio_occ]


        aln_file.close()
        sys.stdout.write("\n")

        compute_results()
