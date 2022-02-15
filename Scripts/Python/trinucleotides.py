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
    for base1 in L_bases:
        for base2 in L_bases:
            for base3 in L_bases:
                if base1 == base2 == base3: # Omit homopolymers
                    continue
                trinucleotide = "-*".join([base1, base2, base3, ""]) # allow gaps within pattern in aligned genome
                pattern += f"({trinucleotide}){{2,}}|"
    pattern = pattern[:-1]
    return pattern


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
        species_name = get_short_name(species_name)
        dictionary_group[species_name] = group_name
        dictionary_color[group_name] = color
    file.close()

    return dictionary_group, dictionary_color


def get_short_name(long_name):
    if len(long_name) <= 10:
        return long_name
    L_name = long_name.replace("_", " ").split(" ")
    L_name[0] = L_name[0][0] + "."
    short_name = " ".join(L_name)
    return short_name

def get_label_name(long_name):
    long_name = long_name.replace("_", " ").split(" ")
    label_name = long_name[0][0] + "."
    for elt in long_name[1:]:
        label_name += " " + elt[:3]
        
    return label_name


def build_group_per_species(filename):
    """
    If no group were defined by user, create N groups (1 for each of the N species)
    """
    dictionary_group = {}
    dictionary_color = {}
    file = open(filename, "r")
    for line in file:
        species_name, color, gc = line.rstrip().split("\t")
        gc = float(gc)
        species_name = get_short_name(species_name)
        group_name = species_name
        dictionary_group[species_name] = group_name
        dictionary_color[group_name] = color
        
    file.close()

    return dictionary_group, dictionary_color



def convert_trinucleotide_to_pattern(sequence: str):
    """
    Converts a trinucleotide to a pattern among: "S-only", "Mostly-S", "W-only", and "Mosly-W"
    S (strong) for bases C or G
    W (weak) for bases A or T
    S and W are relative to number of hydrogen bounds in DNA (3 for GC ; 2 for AT)
    """
    count_CG = sequence.count("C") + sequence.count("G")
    if count_CG == 0:
        pattern = "S only"
    elif count_CG == 1:
        pattern = "Mainly S"
    elif count_CG == 2:
        pattern = "Mainly W"
    else:
        pattern = "W only"
    return pattern



def initiate_accuracy_dict():
    """Initialize dictionary that will store every information needed to plot both
    - distribution of trinucleotides depending on category and length
    - and sequencing accuracy
    For each species group, each trinucleotide type and each length, will store two lists:
        - error rate for each species
        - occurrences for each species
    """
    L_patterns = ["S only", "Mainly S", "Mainly W", "W only"]
    dictionary = {}
    for group_name in L_groups:
        dictionary[group_name] = {}
        for pattern in L_patterns:
            dictionary[group_name][pattern] = {}
            for length in L_lengths:
                dictionary[group_name][pattern][length] = [[], []] # error rate, occurrences
    return dictionary


def initiate_temp_dict():
    """Initialize TEMPORARY dictionary that will store every information needed to plot both
    - distribution of trinucleotides depending on category and length
    - and sequencing accuracy
    For each trinucleotide type and each length, store a list of 3 integers:
        - number of erroneously sequenced bases
        - number of sequenced bases
        - occurrences of such trinucleotide
    """
    L_patterns = ["S only", "Mainly S", "Mainly W", "W only"]
    dictionary = {}
    for pattern in L_patterns:
        dictionary[pattern] = {}
        for length in L_lengths:
            dictionary[pattern][length] = [0, 0, 0]
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
    sw_pattern = convert_trinucleotide_to_pattern(category)
    # Store abundance
    dict_tmp_occurrences_accuracy[sw_pattern][trinucleotide_genome_length][2] += 1
    # Store number of erroneous bases
    nb_errors = 0
    for i, base_genome in enumerate(g):
        base_read = r[i]
        if base_read != base_genome:
            nb_errors += 1
    dict_tmp_occurrences_accuracy[sw_pattern][trinucleotide_genome_length][0] += nb_errors
    dict_tmp_occurrences_accuracy[sw_pattern][trinucleotide_genome_length][1] += len(g)


def compute_results():
    """Compute graph and raw output for current analysis results
    """

    dict_color_trinucl = {"S only": '#8bbeda',
                          "Mainly S": '#b2df8a',
                          "Mainly W": '#fdbf6f',
                          "W only": '#ea484b'}
    patch_S_only = mpatches.Patch(color=dict_color_trinucl["S only"], label='S only')
    patch_main_S = mpatches.Patch(color=dict_color_trinucl["Mainly S"], label='Mainly S')
    patch_main_W = mpatches.Patch(color=dict_color_trinucl["Mainly W"], label='Mainly W')
    patch_W_only = mpatches.Patch(color=dict_color_trinucl["W only"], label='W only')

    nb_groups = len(L_groups)
    nb_lengths = len(L_lengths)



    # --- Abundance of trinucleotides ---
    fig, ax = plt.subplots()
    fig.subplots_adjust(bottom=0.2)

    # just creating another line for x axis labels
    newax = ax.twiny()
    newax.set_frame_on(True)
    newax.patch.set_visible(False)
    newax.xaxis.set_ticks_position('bottom')
    newax.xaxis.set_label_position('bottom')
    newax.spines['bottom'].set_position(('outward', 75))

    fig.set_dpi(300)
    
    # --- If groups have been defined by user: boxplot ---
    if os.path.exists(FILE_SPECIES_GROUP):
        width = 0.5
        offset_trinucleotide = 0.1 + width # space between 2 trinucleotides types for same species group, same length
        offset_group = width # space between results of 2 species groups
        offset_length = 2*width # space between results of 2 different length
        L_start_pos = np.arange(0,
                        (nb_lengths) * (nb_groups * (4*offset_trinucleotide + offset_group) + offset_length),
                        nb_groups * (4*offset_trinucleotide + offset_group) + offset_length)
        for i, group_name in enumerate(dict_occurrences_accuracy):
            for trinucleotide in ["S only", "Mainly S", "Mainly W", "W only"]:
                list_abundance = []
                for length in L_lengths:
                    # -1 values denote that they were absent in dataset
                    # here, remove -1 to keep only "true" values
                    abundances = dict_occurrences_accuracy[group_name][trinucleotide][length][1]
                    abundances = [elt for elt in abundances if elt!=-1]
                    list_abundance += [abundances]
    
                bplot = ax.boxplot(list_abundance, widths=width,
                                   positions=L_start_pos,
                                   patch_artist=True,
                                   zorder=-1)
                colors = [dict_color_trinucl[trinucleotide]]*5
    
                for patch, color in zip(bplot['boxes'], colors):
                    patch.set_facecolor(color)
                    patch.set_edgecolor("black")
                    patch.set_color(color)
                [item.set_color('k') for item in bplot['medians']] # set median line black
    
                L_start_pos += offset_trinucleotide
            L_start_pos += offset_group
    
    
        ax.legend(handles=[patch_S_only, patch_main_S, patch_main_W, patch_W_only],
                  title="Trinucleotides types:", ncol=4)
    
        
    # --- Otherwise: plot ---
    else:
        width = 0.41
        offset_trinucleotide = 0.1 + width # space between 2 trinucleotides types for same species group, same length
        offset_group = width # space between results of 2 species groups
        offset_length = 2*width # space between results of 2 different length
        L_start_pos = np.arange(0,
                        (nb_lengths) * (nb_groups * (4*offset_trinucleotide + offset_group) + offset_length),
                        nb_groups * (4*offset_trinucleotide + offset_group) + offset_length)
        for i, group_name in enumerate(dict_occurrences_accuracy):
            for trinucleotide in ["S only", "Mainly S", "Mainly W", "W only"]:
                list_abundance = []
                for length in L_lengths:
                    abundances = dict_occurrences_accuracy[group_name][trinucleotide][length][1]
                    if abundances == [] or abundances == [-1]:
                        abundances = [float("-inf")]
                    list_abundance += [abundances[0]]
    
                if list_abundance[0] == []:
                    continue
                
                if i == 0:
                    plt.plot(L_start_pos, list_abundance, "o", color=dict_color_trinucl[trinucleotide],
                             label=trinucleotide, ms=4)
                else:
                    plt.plot(L_start_pos, list_abundance, "o", color=dict_color_trinucl[trinucleotide],
                             ms=4)
    
                L_start_pos += offset_trinucleotide
            L_start_pos += offset_group
    
        plt.legend(title="Trinucleotides types:", ncol=4)
    
    #plt.yscale("log")
    plt.xlabel("Number of repetitions")
    ax.set_ylabel("Occurrences ratio of sequenced genomic trinucleotides (%)")

    # First x-axis labels = species groups
    L_xlabels1_pos = []
    p = 0
    for l in L_lengths:
        for g in L_groups:
            L_xlabels1_pos += [p + 1.5 * offset_trinucleotide] # 1.5 to be centered between the 4 categories (=> 3 gaps / 2)
            p += 4 * offset_trinucleotide + offset_group
        p += offset_length
    if os.path.exists(FILE_SPECIES_GROUP):
        L_xlabels1_labels = L_groups * nb_lengths
        ax.set_xticks(L_xlabels1_pos)
        ax.set_xticklabels(L_xlabels1_labels)
    else:
        L_xlabels1_labels = [get_label_name(elt) for elt in L_groups] * nb_lengths 
        ax.set_xticks(L_xlabels1_pos)
        ax.set_xticklabels(L_xlabels1_labels, rotation=90, fontsize=8)
    min_xlim, max_xlim = ax.get_xlim()
    ax.set_xlim(min_xlim-1, max_xlim+2)
    min_ylim, max_ylim = ax.get_ylim()
    ax.set_ylim(max(-0.5, min_ylim-5),
                min(101, max_ylim+5))


    # Second x-axis labels = trinucleotide lengths
    newax.set_xlim(ax.get_xlim()) # set same limits than first x-axis
    L_xlabels2_pos = [np.mean([L_xlabels1_pos[i:i+nb_groups]]) for i in range(0, len(L_xlabels1_pos), nb_groups)]
    L_xlabels2_labels = L_lengths
    newax.set_xticks(L_xlabels2_pos)
    newax.set_xticklabels(L_xlabels2_labels)
    
    plt.savefig(OUTPUT_PLOT + "trinucleotide_occurrences.png")
    plt.close()


    # --- Sequencing accuracy of trinucleotides ---
    # exactly same plot than previously, taking first values in dict instead of second ones

    L_start_pos = np.arange(0,
                            (nb_lengths) * (nb_groups * (4*offset_trinucleotide + offset_group) + offset_length),
                            nb_groups * (4*offset_trinucleotide + offset_group) + offset_length)


    fig, ax = plt.subplots()
    fig.subplots_adjust(bottom=0.2)

    # just creating another line for x axis labels
    newax = ax.twiny()
    newax.set_frame_on(True)
    newax.patch.set_visible(False)
    newax.xaxis.set_ticks_position('bottom')
    newax.xaxis.set_label_position('bottom')
    newax.spines['bottom'].set_position(('outward', 75))

    fig.set_dpi(300)
    # --- If groups have been defined by user: boxplot ---
    if os.path.exists(FILE_SPECIES_GROUP):
        for i, group_name in enumerate(dict_occurrences_accuracy):
            for trinucleotide in ["S only", "Mainly S", "Mainly W", "W only"]:
                list_accuracy = []
                for length in L_lengths:
                    # -1 values denote that they were absent in dataset
                    # here, remove -1 to keep only "true" values
                    accuracy = dict_occurrences_accuracy[group_name][trinucleotide][length][0]
                    accuracy = [elt for elt in accuracy if elt!=-1]
                    list_accuracy += [accuracy]
    
                bplot = ax.boxplot(list_accuracy, widths=width,
                                   positions=L_start_pos,
                                   patch_artist=True,
                                   zorder=-1)
                colors = [dict_color_trinucl[trinucleotide]]*5
    
                for patch, color in zip(bplot['boxes'], colors):
                    patch.set_facecolor(color)
                    patch.set_edgecolor("black")
                    patch.set_color(color)
                [item.set_color('k') for item in bplot['medians']] # set median line black
    
                L_start_pos += offset_trinucleotide
            L_start_pos += offset_group
    
        ax.legend(handles=[patch_S_only, patch_main_S, patch_main_W, patch_W_only],
                  title="Trinucleotides types:", ncol=8)
    
    
        # First x-axis labels = species groups
        L_xlabels1_pos = []
        p = 0
        for l in L_lengths:
            for g in L_groups:
                L_xlabels1_pos += [p + 1.5 * offset_trinucleotide] # 1.5 to be centered between the 4 categories (=> 3 gaps / 2)
                p += 4 * offset_trinucleotide + offset_group
            p += offset_length
        L_xlabels1_labels = L_groups * nb_lengths
        ax.set_xticks(L_xlabels1_pos)
        ax.set_xticklabels(L_xlabels1_labels)
    
    
        # Second x-axis labels = trinucleotide lengths
        newax.set_xlim(ax.get_xlim()) # set same limits than first x-axis
        L_xlabels2_pos = [np.mean([L_xlabels1_pos[i:i+nb_groups]]) for i in range(0, len(L_xlabels1_pos), nb_groups)]
        L_xlabels2_labels = L_lengths
        newax.set_xticks(L_xlabels2_pos)
        newax.set_xticklabels(L_xlabels2_labels)
    
        
    # --- Otherwise: plot ---
    else:
        for i, group_name in enumerate(dict_occurrences_accuracy):
            for trinucleotide in ["S only", "Mainly S", "Mainly W", "W only"]:
                list_accuracy = []
                for length in L_lengths:
                    accuracy = dict_occurrences_accuracy[group_name][trinucleotide][length][0]
                    if accuracy == [] or accuracy == [-1]:
                        accuracy = [float("-inf")]
                    list_accuracy += [accuracy[0]]

                if list_accuracy[0] == []:
                    continue
                
                if i == 0:
                    plt.plot(L_start_pos, list_accuracy, "o", color=dict_color_trinucl[trinucleotide],
                             label=trinucleotide, ms=4)
                else:
                    plt.plot(L_start_pos, list_accuracy, "o", color=dict_color_trinucl[trinucleotide],
                             ms=4)
    
                L_start_pos += offset_trinucleotide
            L_start_pos += offset_group
    
        plt.legend(title="Trinucleotides types:", ncol=4)
    
    
    # First x-axis labels = species groups
    L_xlabels1_pos = []
    p = 0
    for l in L_lengths:
        for g in L_groups:
            L_xlabels1_pos += [p + 1.5 * offset_trinucleotide] # 1.5 to be centered between the 4 categories (=> 3 gaps / 2)
            p += 4 * offset_trinucleotide + offset_group
        p += offset_length
    
    if os.path.exists(FILE_SPECIES_GROUP):
        L_xlabels1_labels = L_groups * nb_lengths
        ax.set_xticks(L_xlabels1_pos)
        ax.set_xticklabels(L_xlabels1_labels)
    else:
        L_xlabels1_labels = [get_label_name(elt) for elt in L_groups] * nb_lengths
        ax.set_xticks(L_xlabels1_pos)
        ax.set_xticklabels(L_xlabels1_labels, rotation=90, fontsize=8)
    min_xlim, max_xlim = ax.get_xlim()
    ax.set_xlim(min_xlim-1, max_xlim+2)
    min_ylim, max_ylim = ax.get_ylim()
    ax.set_ylim(max(-0.5, min_ylim-5),
                min(101, max_ylim+5))

    # Second x-axis labels = trinucleotide lengths
    newax.set_xlim(ax.get_xlim()) # set same limits than first x-axis
    L_xlabels2_pos = [np.mean([L_xlabels1_pos[i:i+nb_groups]]) for i in range(0, len(L_xlabels1_pos), nb_groups)]
    L_xlabels2_labels = L_lengths
    newax.set_xticks(L_xlabels2_pos)
    newax.set_xticklabels(L_xlabels2_labels)
    
    plt.xlabel("Number of repetitions")
    ax.set_ylabel("Correctly sequenced trinucleotides (%)")
    plt.savefig(OUTPUT_PLOT + "trinucleotide_accuracy.png")
    plt.close()


    # --- Saving results in file ---
    dict_trinucl_occ = {}
    dict_trinucl_accuracy = {}
    for group_name in dict_occurrences_accuracy:
        dict_trinucl_occ[group_name] = {}
        dict_trinucl_accuracy[group_name] = {}
        for trinucleotide_category in dict_occurrences_accuracy[group_name]:
            dict_trinucl_occ[group_name][trinucleotide_category] = {}
            dict_trinucl_accuracy[group_name][trinucleotide_category] = {}
            for length in dict_occurrences_accuracy[group_name][trinucleotide_category]:
                list_accuracy, list_occ = dict_occurrences_accuracy[group_name][trinucleotide_category][length]
                dict_trinucl_occ[group_name][trinucleotide_category][length] = list_occ
                dict_trinucl_accuracy[group_name][trinucleotide_category][length] = list_accuracy


    output = open(OUTPUT_RAW + "trinucleotides_occurrence_accuracy.txt", "w")
    write_dictionary(dict_trinucl_occ, "Occurences", output)
    write_dictionary(dict_trinucl_accuracy, "\nAccuracy", output)
    output.close()


def write_dictionary(dictionary, text, file):
    """
    Write dictionaries to output file
    """
    file.write(text + "\n")
    for group_name in dictionary:
        file.write(group_name + "\n")
        for category in dictionary[group_name]:
            file.write(category + "\n")
            for length in dictionary[group_name][category]:
                results = [elt if elt!=-1 else "NaN" for elt in dictionary[group_name][category][length]]
                results = "\t".join(map(str, results))
                file.write("\t" + f"{length} {results}" +  "\n")
                #file.write("\t" + length, dictionary[group_name][category][length], "\n")
        file.write("\n")




# --------------------------------------------------------------------------------------------------
# Main

if __name__ == "__main__":


    ## Parses arguments
    NUMBER_EXPECTED_ARGUMENTS = 9
    if len(sys.argv) != NUMBER_EXPECTED_ARGUMENTS + 1:
        print(f"ERROR: Wrong number of arguments: {NUMBER_EXPECTED_ARGUMENTS} expected but {len(sys.argv)-1} given.")
        sys.exit(2)
    ALN_EXPL_DIRNAME, OUTPUT_RAW, OUTPUT_PLOT, FILE_SPECIES_GC, FILE_SPECIES_GROUP, MIN_LENGTH, MAX_LENGTH, THRESHOLD_OCCURENCES, NB_MAX_ALN = sys.argv[1:]

    MIN_LENGTH = int(MIN_LENGTH)
    MAX_LENGTH = int(MAX_LENGTH)
    THRESHOLD_OCCURENCES = int(THRESHOLD_OCCURENCES) # Minimal number of occurrences to take into account the result
    NB_MAX_ALN = int(NB_MAX_ALN) # Maximal number of alignment processed, to speed up

    if NB_MAX_ALN == -1:
        NB_MAX_ALN = float("inf")

    # --- Parameters ---

    # Some global variables
    L_bases = ["A", "C", "G", "T"]
    pattern = get_regex_pattern() # regex pattern for searching for trinucleotides

    # Trinucleotide "length" defined as their number of bases
    #   for example ACGACG is of length 6 (thus can only be multiples of 3)
    L_lengths = list(np.arange(MIN_LENGTH, MAX_LENGTH, 3)) + [f"{MAX_LENGTH}+"]


    # Group species:
    #   - by user defined group (if exists)
    if os.path.exists(FILE_SPECIES_GROUP):
        dict_species_group, dict_species_group_color = get_groups(FILE_SPECIES_GROUP)
    #   - else build artificial groups: one per species
    else:
        dict_species_group, dict_species_group_color = build_group_per_species(FILE_SPECIES_GC)
    L_groups = sorted(set(dict_species_group.values()))
    
    
    # Initiate dictionary to store results
    dict_occurrences_accuracy = initiate_accuracy_dict()

    # --- Analyse alignment files ---
    for aln_filename in os.listdir(ALN_EXPL_DIRNAME):

        aln_path = ALN_EXPL_DIRNAME + aln_filename
        NB_TOT_ALN = get_total_number_of_lines(aln_path) / 3
        STARTING_TIME = time.time()

        species_name = aln_filename.split(".txt")[0]
        print("  ", species_name)
        group_name = dict_species_group[get_short_name(species_name)]
        aln_file = open(aln_path, "r")

        # temporary dictionary storing result for current species
        #    dict[SW-pattern][length] = [nb_err, nb_seq, nb_occ] ;  with :
        #       - nb_err = number of erroneously sequenced bases
        #       - nb_seq = number of sequenced bases
        #       - nb_occ = occurrences of such trinucleotide
        dict_tmp_occurrences_accuracy = initiate_temp_dict()

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

            if nb_aln_done >= NB_MAX_ALN:
                break


            # --- Compute sequencing errors associated with genomic trinucleotides ---
            start_all, end_all = -1, -1
            for res in re.finditer(pattern, genome_aln):
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


                trinucleotide_genome_length = len(trinucleotide_genome.replace("-", ""))
                if trinucleotide_genome_length >= MAX_LENGTH:
                    trinucleotide_genome_length = f"{MAX_LENGTH}+"
                start_read, end_read = try_extend_read_trinucleotide(start, end, read_aln)
                start_all = min(start_genome, start_read)
                end_all = max(end_genome, end_read)

                # Update dictionary storing results
                update_results(genome_aln[start_all: end_all],
                               read_aln[start_all: end_all])


        # Merge temporary result to final dictionary storing results
        for sw_pattern in dict_tmp_occurrences_accuracy:
            for length in dict_tmp_occurrences_accuracy[sw_pattern]:
                nb_err, nb_tot, nb_occ = dict_tmp_occurrences_accuracy[sw_pattern][length]
                if nb_occ > THRESHOLD_OCCURENCES:
                    ratio_occ = nb_occ / total_nb_sequenced_bases * 100
                    accuracy_rate = round(100 - (nb_err / nb_tot * 100), 2)
                    dict_occurrences_accuracy[group_name][sw_pattern][length][0] += [accuracy_rate]
                    dict_occurrences_accuracy[group_name][sw_pattern][length][1] += [ratio_occ]
                else:
                    dict_occurrences_accuracy[group_name][sw_pattern][length][0] += [-1]
                    dict_occurrences_accuracy[group_name][sw_pattern][length][1] += [-1]



        aln_file.close()
        sys.stdout.write("\n")


        compute_results()


