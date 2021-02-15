#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Analyses substitution errors and output results as a plot
Results are store separating low GC bacteria, high GC bacteria and human data
"""

# --------------------------------------------------------------------------------------------------
# Packages

import sys
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

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
# Main

if __name__ == "__main__":

    ## Parses arguments
    if len(sys.argv) != 4:
        print(f"ERROR: Wrong number of arguments: 3 expected but {len(sys.argv)-1} given.")
        sys.exit(2)
    ERRORS_DIRNAME = sys.argv[1]
    OUTPUT_GRAPH = sys.argv[2]
    FILENAME_SPECIES_GC_COLOR = sys.argv[3]

    if ERRORS_DIRNAME[-1] != "/":
        ERRORS_DIRNAME += "/"
    if OUTPUT_GRAPH[-1] != "/":
        OUTPUT_GRAPH += "/"



    # Initialize dictionary that will store sustitution error rates
    dict_substitutions = {}
    for base1 in ["A", "C", "T", "G"]:
        for base2 in ["A", "C", "T", "G"]:
            if base1 == base2:
                continue
            substitution = base1 + base2
            dict_substitutions[substitution] = {} # will store occurrences for each species

    # Browse results for datasets
    print(f"Browsing files in {ERRORS_DIRNAME}")
    for filename in os.listdir(ERRORS_DIRNAME):

        species_name = filename.split(".")[0]

        set_species_name = set()
        for substitution in dict_substitutions:
            if species_name not in dict_substitutions[substitution]:
                dict_substitutions[substitution][species_name] = 0
            set_species_name.add(species_name)

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
                    if substitution in dict_substitutions:
                        dict_substitutions[substitution][species_name] += occurrences
                        sum_substitution += occurrences
                break
        file.close()

        # convert substitution occurrences to ratio (compared to all other for this species)
        for substitution in dict_substitutions:
            occurrences = dict_substitutions[substitution][species_name]
            ratio = round(occurrences / sum_substitution * 100, 2)
            dict_substitutions[substitution][species_name] = ratio


    # Gather results in categories: low/high bacteria and human
    print("Gathering data")
    ## get gc information for each species
    dict_species_gc = {}
    species_gc_file = open(FILENAME_SPECIES_GC_COLOR, "r")
    for line in species_gc_file:
        species_name, _, gc = line.rstrip().split(" ; ")
        species_name = species_name.replace(" ", "_")
        if "Human" in species_name:
            dict_species_gc[species_name] = "human"
            continue
        gc = float(gc)
        if gc > 50:
            dict_species_gc[species_name] = "high GC bacteria"
        else:
            dict_species_gc[species_name] = "low GC bacteria"
    species_gc_file.close()

    ## gather data
    dict_to_plot = {"low GC bacteria": [], "high GC bacteria": [], "human": []}
    list_ordered_substi = ["AC", "CA", "AG", "GA", "AT", "TA", "CG", "GC", "CT", "TC", "GT", "TG"]
    for substitution in list_ordered_substi:
        list_temp_low, list_temp_high, list_temp_human = [], [], []
        for species_name in dict_substitutions[substitution]:
            gc_category = dict_species_gc[species_name]
            substitution_rate = dict_substitutions[substitution][species_name]
            if "low" in gc_category:
                list_temp_low += [substitution_rate]
            elif "high" in gc_category:
                list_temp_high += [substitution_rate]
            else:
                list_temp_human += [substitution_rate]
        dict_to_plot["low GC bacteria"] += [list_temp_low]
        dict_to_plot["high GC bacteria"] += [list_temp_high]
        dict_to_plot["human"] += [list_temp_human]


    # Plot
    print("Constructing the plot")
    fig = plt.figure()
    fig.set_dpi(300.0)
    boxplot_width = 0.2
    list_labels = []
    for substitution in list_ordered_substi:
        list_labels += [substitution[0] + "-" + substitution[1]]
    list_pos_low, list_pos_high, list_pos_human = [], [], []
    diff_inter = 1
    diff_intra = 0.8
    diff_gc_category = 0.22
    pos = -diff_inter
    for i in range(6):
        pos += diff_inter
        list_pos_low += [pos - diff_gc_category]
        list_pos_high += [pos]
        list_pos_human += [pos + diff_gc_category]
        pos += diff_intra
        list_pos_low += [pos - diff_gc_category]
        list_pos_high += [pos]
        list_pos_human += [pos + diff_gc_category]

    box_low = plt.boxplot(dict_to_plot["low GC bacteria"], positions=list_pos_low, labels=[""]*len(list_pos_low),
                patch_artist=True, widths=boxplot_width)
    color_low = "#4477AA"
    for item in ['boxes', 'whiskers', 'fliers', 'caps']:
        plt.setp(box_low[item], color=color_low)
    for item in ['medians']:
        plt.setp(box_low[item], color="black")

    box_high = plt.boxplot(dict_to_plot["high GC bacteria"], positions=list_pos_high, labels=list_labels,
                patch_artist=True, widths=boxplot_width)
    color_high = "#EE6677"
    for item in ['boxes', 'whiskers', 'fliers', 'caps']:
        plt.setp(box_high[item], color=color_high)
    for item in ['medians']:
        plt.setp(box_high[item], color="black")

    color_human = "#CCBB44"
    for i, list_res_substi in enumerate(dict_to_plot["human"]):
        for res_substi in list_res_substi:
            plt.plot(list_pos_human[i], res_substi, "x", ms=12, color=color_human)

    plt.xticks(ticks=list_pos_high, labels=list_labels)
    plt.xlim(-0.5, 10.3)

    plt.xlabel("Substitution (genomic base - read base)")
    plt.ylabel("Frequency (%)")

    legend_colors = [color_low, color_high, color_human]
    legend_texts = ["Low GC bacteria", "High GC bacteria", "Human"]
    legend_patches = [ mpatches.Patch(color=legend_colors[i], label="{:s}".format(legend_texts[i])) \
                      for i in range(len(legend_texts)) ]
    plt.legend(handles=legend_patches, title="Species:", ncol=3)

    plt.savefig(OUTPUT_GRAPH + "substitution_errors.png")

    plt.close()
