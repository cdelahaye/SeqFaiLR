#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Analyses substitution errors and output results as a plot
Results are displayed grouping reads according to their group, as defined by the user during
  initialization step
Or detailled for each species
"""

# --------------------------------------------------------------------------------------------------
# Packages

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

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

def init_dict_substitutions():
    """
    Initialize a dictionary that will store substitution error occurrences for each species.
    Dictionary of form: dictionary[substitution_type][species_name] = substitution_occurrences

    """
    L_bases = ["A", "C", "G", "T"]
    L_substitutions = [base1+base2 for base1 in L_bases for base2 in L_bases if base1!=base2]
    dictionary = {}
    for substitution in L_substitutions:
        dictionary[substitution] = {} # will store occurrences for each species
    return dictionary

def get_short_name(long_name):
    L_name = long_name.replace("_", " ").split(" ")
    L_name[0] = L_name[0][0] + "."
    short_name = " ".join(L_name)
    return short_name


def get_substitution_errors(dictionary, dirname):
    """
    Browse errors for each species (information in files of dirname)
    Store results in dictionary
    """
    for filename in os.listdir(ERRORS_DIRNAME):

        species_name = filename.split(".")[0]

        for substitution in dictionary:
            if species_name not in dictionary[substitution]:
                dictionary[substitution][species_name] = 0

        # extract information from files, for each species
        file = open(ERRORS_DIRNAME + filename, "r")
        sum_substitution = 0
        for line in file:
            if "Substitutions" in line and "Occurrences" in line:
                while True:
                    line = file.readline().rstrip()
                    if not line:
                        break
                    substitution, occurrences = line.split("\t")
                    occurrences = int(occurrences)
                    if substitution in dictionary:
                        dictionary[substitution][species_name] += occurrences
                        sum_substitution += occurrences
                break
        file.close()

        # convert substitution occurrences to ratio (compared to all others for this species)
        for substitution in dictionary:
            occurrences = dictionary[substitution][species_name]
            ratio = round(occurrences / sum_substitution * 100, 2)
            dictionary[substitution][species_name] = ratio

    return dictionary


def get_species_groups(filename_groups, filename_gc):
    """
    Parse filename_groups to get user-defined species groups.
    If not defined, create a group for each species
    Outputs a dictionary grouping species.
    """
    dictionary = {}

    # get user defined groups
    if os.path.exists(filename_groups):
        file_groups = open(filename_groups, "r")
        while True:
            line = file_groups.readline()
            if not line:
                break
            species_name = line.split("\t")[1].split(".")[0]
            group_name = line.split("\t")[0]
            dictionary[species_name] = group_name
        file_groups.close()

    # otherwise, define a group for each species
    else:
        file_gc = open(filename_gc, "r")
        while True:
            line = file_gc.readline()
            if not line:
                break
            species_name, _, _ = line.split()
            dictionary[species_name] = get_short_name(species_name)
        file_gc.close()

    return dictionary


def gather_data_in_groups(dictionary_groups, dictionary_substitutions):
    """
    Gather substitution results in each group
    Output in dictionary that will be used to plot final results
    """
    dictionary = {}

    # ordered substitution for graph output:
    L_ordered_substitutions = ["AC", "CA", "AG", "GA", "AT", "TA",
                               "CG", "GC", "CT", "TC", "GT", "TG"]

    # gather results for each group for each substitution
    for substitution in L_ordered_substitutions:
        dict_tmp = {}
        for species_name in dictionary_substitutions[substitution]:
            group_name = dictionary_groups[species_name]
            substitution_rate = dictionary_substitutions[substitution][species_name]
            if group_name not in dictionary:
                dictionary[group_name] = []
            if group_name not in dict_tmp:
                dict_tmp[group_name] = []
            dict_tmp[group_name] += [substitution_rate]
        for group_name in dict_tmp:
            dictionary[group_name] += [dict_tmp[group_name]]

    return dictionary


def plot_substitution_errors(dictionary, dictionary_groups):
    """
    Plot results contained in dictionary
    """

    # Substitution labels
    L_ordered_substitutions = ["AC", "CA", "AG", "GA", "AT", "TA",
                               "CG", "GC", "CT", "TC", "GT", "TG"]
    substitution_number = int(len(L_ordered_substitutions) / 2)
    list_labels = []
    for substitution in L_ordered_substitutions:
        list_labels += [substitution[0] + "-" + substitution[1]]

    # --- Plot parameters ---
    # general plot
    fig = plt.figure()
    fig.set_dpi(500.0)
    # spacing between boxplots
    group_number = len(dictionary)
    offset_inter_substitutions = 3 # offset between 2 substitution groups (e.g. AC/CA and AG/GA)
    offset_intra_substitutions = 2 # offset between 2 substitutions (e.g. AC and CA)
    offset_group = 1 # offset between species categories within a substitution result
    position = 0 # initiate position
    dict_positions = {} # dictionary to store boxplot positions
    # group colors
    list_colors = sns.diverging_palette(240, 10, n=group_number, as_cmap=False)
    list_colors_hex = list_colors.as_hex() # get hexadecimal code for each color

    # Get group names and create dictionary to store their boxplot positions
    set_group_names = set(dictionary_groups.values())
    for group_name in set_group_names:
        dict_positions[group_name] = []

    # Compute boxplot positions and label positions
    L_label_positions = []
    for i in range(substitution_number):
        tmp_pos = []
        for group_name in set_group_names:
            dict_positions[group_name] += [position]
            tmp_pos += [position]
            position += offset_group
        L_label_positions += [np.median(tmp_pos)]
        position -= offset_group
        position += offset_intra_substitutions
        tmp_pos = []
        for group_name in set_group_names:
            dict_positions[group_name] += [position]
            tmp_pos += [position]
            position += offset_group
        L_label_positions += [np.median(tmp_pos)]
        position -= offset_group
        position += offset_inter_substitutions
    
    # --- If user defineds groups: boxplots ---
    if os.path.exists(FILENAME_SPECIES_GROUPS):
        for i, group_name in enumerate(set_group_names):
            tmp_boxplot = plt.boxplot(dictionary[group_name], positions=dict_positions[group_name],
                                      labels=[""]*len(dict_positions[group_name]), patch_artist=True)
            tmp_color = list_colors_hex[i]
            for item in ["boxes", "whiskers", "fliers", "caps"]:
                plt.setp(tmp_boxplot[item], color=tmp_color)
            plt.setp(tmp_boxplot["medians"], color="black")
    
        legend_colors = list_colors_hex
        legend_texts = list(set_group_names)
        legend_patches = [ mpatches.Patch(color=legend_colors[i], label="{:s}".format(legend_texts[i])) \
                          for i in range(len(legend_texts)) ]
        plt.legend(handles=legend_patches, title="Species groups:", ncol=2)
    
        
    # --- Otherwise, plot for each species ---
    else:
        for i, group_name in enumerate(set_group_names):
            color = list_colors_hex[i]
            plt.plot(dict_positions[group_name], dictionary[group_name], "o", ms=4, 
                     color=color, label=group_name)
        plt.legend(title="Species groups:", ncol=2)

    plt.xticks(ticks=L_label_positions, labels=list_labels)
    plt.xlim(-1, position)

    plt.xlabel("Substitution (genomic base - read base)")
    plt.ylabel("Frequency (%)")

    plt.savefig(OUTPUT_GRAPH)

    plt.close()




# --------------------------------------------------------------------------------------------------
# Main

if __name__ == "__main__":

    # Parses arguments
    NUMBER_EXPECTED_ARGUMENTS = 4
    if len(sys.argv) != NUMBER_EXPECTED_ARGUMENTS + 1:
        print(f"ERROR: Wrong number of arguments: {NUMBER_EXPECTED_ARGUMENTS} expected but {len(sys.argv)-1} given.")
        sys.exit(2)
    ERRORS_DIRNAME, OUTPUT_GRAPH, FILENAME_SPECIES_GC_COLOR, FILENAME_SPECIES_GROUPS = sys.argv[1:]

    if ERRORS_DIRNAME[-1] != "/":
        ERRORS_DIRNAME += "/"


    # Initialize dictionary that will store sustitution error rates
    dict_substitutions = init_dict_substitutions()


    # Browse results for each species
    dict_substitutions = get_substitution_errors(dict_substitutions, ERRORS_DIRNAME)


    # Gather results in groups, if user defined them
    dict_species_group = get_species_groups(FILENAME_SPECIES_GROUPS, FILENAME_SPECIES_GC_COLOR,)
    dict_to_plot = gather_data_in_groups(dict_species_group, dict_substitutions)


    # Plot results
    plot_substitution_errors(dict_to_plot, dict_species_group)

