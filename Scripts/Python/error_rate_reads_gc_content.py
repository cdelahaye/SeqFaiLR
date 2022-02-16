#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
For each group, computes mean error rates (total and detailled)
depending on the GC content of the read
"""


# --------------------------------------------------------------------------------------------------
# Packages

import sys
import os
import matplotlib.pyplot as plt

params = {'legend.fontsize': 16,
          'legend.title_fontsize': 16,
          'legend.labelspacing': 0.1,
          'legend.borderpad': 0.3,
          'legend.columnspacing': 1,
          'legend.handletextpad': 0.5,
          'legend.handlelength': 1.2,
          'figure.figsize': (12, 8),
          'axes.labelsize': 16,
          'axes.titlesize': 16,
          'xtick.labelsize': 16,
          'ytick.labelsize': 16}
plt.rcParams.update(params)

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

def get_gc_color(filename):
    """
    Returns a dictionary of species names (keys), associated with a color (value) attributed
    according to its GC content
    Also compute a list of increasing GC content species names
    """
    dictionary = {}
    dict_gc_species = {}
    file = open(filename, "r")
    for line in file:
        species_name, color, gc = line.rstrip().split("\t")
        dictionary[species_name] = color
        if gc not in dict_gc_species:
            dict_gc_species[gc] = []
        dict_gc_species[gc] += [species_name]
    file.close()

    list_ordered_species_gc = []
    for gc in sorted(dict_gc_species):
        list_ordered_species_gc += [elt for elt in dict_gc_species[gc]]

    return dictionary, list_ordered_species_gc





def get_error_rates(genome, read) :
    """ Computes error rates (total, mismatch, insertion and deletion) for a given alignment
    between a genome and and read
    """
    mismatches, insertions, deletions = 0, 0, 0
    length = len(genome.replace("N", ""))
    for i, base_genome in enumerate(genome) :
        base_read = read[i]
        if base_genome == "N" :
            continue
        elif base_genome == "-" :
            insertions += 1
        elif base_read == "-" :
            deletions += 1
        elif base_genome != base_read :
            mismatches += 1
    total = round((mismatches + insertions + deletions) / length * 100, 1)
    mismatches = round(mismatches / length * 100, 4)
    insertions = round(insertions / length * 100, 4)
    deletions = round(deletions / length * 100, 4)
    return [total, mismatches, insertions, deletions]


def get_gc_content_of_read(read):
    """ Computes the GC content of a given read
    """
    r = read.replace("-", "")
    gc_content = round((r.count("G") + r.count("C")) / len(r) * 100)
    return gc_content



def get_error_rates_gc():
    """
    Parse alignment files to get error rates depending on GC content of reads
    Store results in dictionary
    """
    dictionary = {}

    for aln_filename in os.listdir(ALN_EXPL_DIRNAME):
        if os.path.isdir(ALN_EXPL_DIRNAME + aln_filename):
            continue

        # Get species name
        species_name = aln_filename.split(".txt")[0]
        print("  ", species_name)

        # --- Parse alignment file ---
        file = open(ALN_EXPL_DIRNAME + aln_filename, "r")

        # initialize dictionary
        dictionary[species_name] = {}
        for error_type in L_error_types:
            dictionary[species_name][error_type] = {}

        aln_count = 0 # count number of alignments processed
        while True:
            header = file.readline()
            if not header:
                break
            genome = file.readline().rstrip()
            read = file.readline().rstrip()

            # Threshold if you want to limit analyses and speed up computation
            aln_count += 1
            if aln_count > NB_MAX_ALN:
                break

            # Get GC content and error rate of read
            gc = get_gc_content_of_read(read)
            L_error_rates = get_error_rates(genome, read)

            # Update dictionary
            for i, error_rate in enumerate(L_error_rates):
                error_type = L_error_types[i]
                if gc not in dictionary[species_name][error_type]:
                    dictionary[species_name][error_type][gc] = [0, 0] # error rate value, occurrences
                dictionary[species_name][error_type][gc][0] += error_rate
                dictionary[species_name][error_type][gc][1] += 1
        file.close()

    return dictionary



def gather_results_groups(dictionary_results_species, dictionary_groups):
    """
    From a dictionary_results_species that stored results for each species, gather these results
    in groups of species, as defined by user
    """
    dictionary = {}
    for species_name in dictionary_results_species:
        group_name = dictionary_groups[species_name]
        if group_name not in dictionary:
            dictionary[group_name] = {}
            for error_type in L_error_types:
                dictionary[group_name][error_type] = {}
        for error_type in L_error_types:
            for gc in dictionary_results_species[species_name][error_type]:
                sum_errors, nb_occ_errors = dictionary_results_species[species_name][error_type][gc]
                if gc not in dictionary[group_name][error_type]:
                    dictionary[group_name][error_type][gc] = [0, 0]
                dictionary[group_name][error_type][gc][0] += sum_errors
                dictionary[group_name][error_type][gc][1] += nb_occ_errors
    return dictionary

def get_short_name(long_name):
    L_name = long_name.replace("_", " ").split(" ")
    if len(L_name) == 1:
        return long name
    L_name[0] = L_name[0][0] + "."
    short_name = " ".join(L_name)
    return short_name

def write_results(dictionary, output_filename):
    """
    Write results in output_filename, as follow:
        - name of species (or group)
        - list (tab separated) of gc content
        - list (tab separated) of associated mean XXX error rate (for each XXX eror type)

    """
    file = open(OUTPUT_RAW + output_filename, "w")
    for entity in dictionary: # can be a species or a group
        file.write(f"{entity}\n")
        for error_type in dictionary[entity]:
            L_gc = []
            L_mean_error_rate = []
            for gc in sorted(dictionary[entity][error_type]):
                sum_errors, nb_occurences = dictionary[entity][error_type][gc]
                L_gc += [str(gc)]
                L_mean_error_rate += [str(round(sum_errors / nb_occurences, 4))]
            file.write("GC\t" + "\t".join(L_gc) + "\n")
            file.write(f"{error_type}\t" + "\t".join(L_mean_error_rate) + "\n")
        file.write("\n")
    file.close()


def plot_results(dictionary, output_filename, dict_colors, L_order, category):
    """
    Plot results contained in dictionary, in output_filename
    """

    for error_type in L_error_types:
        # initialize plot
        fig = plt.figure()
        ax = fig.add_subplot(111)
        fig.set_dpi(300.0)
        box = ax.get_position()
        ax.set_position([0.1, 0.12, box.width*0.8, box.height])

        # add results for each species
        for entity in dictionary:
            color = dict_colors[entity]
            L_gc, L_error_rate = [], []
            for gc in sorted(dictionary[entity][error_type]):
                sum_errors, nb_occurrences = dictionary[entity][error_type][gc]
                error_rate = sum_errors / nb_occurrences
                L_gc += [gc]
                L_error_rate += [error_rate]
            ax.plot(L_gc, L_error_rate, color=color, label=entity)

        # Labels and legend
        ax.set(xlabel='Average GC content of read (%)',
               ylabel=f"{error_type} error rate (%)")
        #   reorder legend labels
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ordered_label_list = [elt for elt in L_order if elt in by_label]
        ordered_label_values = [by_label[k] for k in ordered_label_list]
        if category == "Species":
            ordered_label_list_short_name = [get_short_name(elt) for elt in ordered_label_list]
        else:
            ordered_label_list_short_name = [elt for elt in ordered_label_list]
        ax.legend(ordered_label_values, ordered_label_list_short_name, title=category,
                  fontsize=12, title_fontsize=14,
                   bbox_to_anchor=(1, 1)) # place legend outside plot
        plt.savefig(OUTPUT_PLOT + output_filename.replace(".png", f"_{error_type}.png"))
        plt.close()




# --------------------------------------------------------------------------------------------------
# Main

if __name__ == "__main__":


    ## Parses arguments
    NUMBER_EXPECTED_ARGUMENTS = 7
    if len(sys.argv) != NUMBER_EXPECTED_ARGUMENTS + 1:
        print(f"ERROR: Wrong number of arguments: {NUMBER_EXPECTED_ARGUMENTS} expected but {len(sys.argv)-1} given.")
        sys.exit(2)
    ALN_EXPL_DIRNAME, OUTPUT_RAW, OUTPUT_PLOT, FILE_SPECIES_GC, FILE_SPECIES_GROUPS, NB_MAX_ALN, MIN_OCC_TO_PLOT= sys.argv[1:]
    NB_MAX_ALN = int(NB_MAX_ALN)
    MIN_OCC_TO_PLOT = int(MIN_OCC_TO_PLOT) # Minimal number of occurences for a certain GC content value to be plot

    if NB_MAX_ALN == -1:
        NB_MAX_ALN = float("inf")


    # List of error categories, and associated colors for plot
    L_error_types = ["Total", "Mismatch", "Insertion", "Deletion"]
    
    # Parse alignments to get error rates depending on GC content
    #   the dictionary will store, for each species, the error rate depending on GC content of reads
    dict_gc_error = get_error_rates_gc()

    # If groups have been defined, compute for groups
    if os.path.exists(FILE_SPECIES_GROUPS):
        # get group of each species
        dict_species_group, dict_species_group_color = get_groups(FILE_SPECIES_GROUPS)

        # parse alignments to get error rates depending on GC content
        #   the dictionary will store, for each species, the error rate depending on GC content of reads
        dict_gc_error_group = gather_results_groups(dict_gc_error, dict_species_group)

        # write results in raw txt file
        write_results(dict_gc_error_group, "error_rate_gc_groups.txt") # for each group

        # plot results
        plot_results(dict_gc_error_group, "error_rate_gc_groups.png", # for each group, ordered alphabetically
                     dict_species_group_color, sorted(set(dict_species_group.values())), "Groups")


    # Otherwise, compute details for each species
    else:
        dict_species_gc_color, L_ordered_species = get_gc_color(FILE_SPECIES_GC)


        # write results in raw txt file
        write_results(dict_gc_error, "error_rate_gc_species.txt") # for each species

        # plot results
        plot_results(dict_gc_error, "error_rate_gc_species.png", # for each species, ordered according GC content
                     dict_species_gc_color, L_ordered_species, "Species")

