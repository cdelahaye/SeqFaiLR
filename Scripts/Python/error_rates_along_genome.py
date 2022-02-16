#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Computes error rates (mismatches, insertions, deletions, and global error rate)
depending on the relative position in the reference genome.
"""

# --------------------------------------------------------------------------------------------------
# Packages

import re
import sys
import glob
import os
import time
import datetime
import math
import numpy as np
import matplotlib.pyplot as plt



# --------------------------------------------------------------------------------------------------
# Global parameters

# Colors for plot:
COLOR_MISMATCH = "#DDAA33"
COLOR_INSERTION = "#BB5566"
COLOR_DELETION = "#004488"
COLOR_AVERAGE = "k"
L_COLORS = [COLOR_MISMATCH, COLOR_INSERTION, COLOR_DELETION, COLOR_AVERAGE]

# Error types
L_errors_types = ["Mismatch", "Insertion", "Deletion", "Total"]


# --------------------------------------------------------------------------------------------------
# Functions

def initialize(L_errors):
    """
    Initialize variables for future plot.

    dictionary: for each type of error (mismatch, insertion, deletion, and all errors),
                save error rates at each relative position in the genome (0 to 99%)
    L_species: list of species names
    Both are empty at the moment
    """
    dictionary = {}
    for error_type in L_errors:
        dictionary[error_type] = {}
        for position in range(100):
            dictionary[error_type][position] = []
    L_species = []

    return dictionary, L_species


def checks_is_forward_alignment(path: str) -> bool:
    """
    Test if given path is an alignment file of forward strand (reverse strand are dealt later)
    Returns a boolean (True/False)
    """
    is_forward = True

    # If the path is not a file -> False
    if not os.path.isfile(path):
        is_forward = False

    # If it is not an alignment file -> False
    elif ".sam" not in path:
        is_forward = False

    # If it is an alignment of reverse strand -> False
    # (they are not ignored but handled later)
    elif "reverse" in path:
        is_forward = False

    # If the reference genome is in multiple chromosomes -> False
    # (this script currently can handle only single chromosome genomes, any suggestions welcome!)
    else:
        species_name = os.path.basename(path).split(".sam")[0]
        print(species_name)
        reference_genome_filename = glob.glob(REF_GEN_DIRNAME + species_name + ".f*")[0]
        count_start_part_fasta = 0
        file = open(reference_genome_filename, "r")
        for line in file:
            if line[0] == ">":
                count_start_part_fasta += 1
            if count_start_part_fasta > 1:
                is_forward = False
                continue

    if not is_forward:
        return False, "", ""
    else:
        return True, species_name, reference_genome_filename


def get_total_number_of_lines(filename: str) -> int:
    """
    Returns number of lines of filename given as input
    """
    with open(filename, "r") as file:
        return sum(1 for _ in file)


def get_reference_genome_length(filename) -> int:
    """
    Computes the total length of reference genome given as argument (expected to be a .fasta file)
    """
    genome_length = 0
    with open(filename, "r") as ref_genome:
        for line in ref_genome:
            if line[0] != ">": # skip
                genome_length += len(line.replace("\n", ""))
    return genome_length


def extract_errors_from_sam(sam_file, nb_tot_alignments, dict_tmp):
    """
    Parse SAM file and extract errors from CIGAR information
    """

    # counters for progressing bar
    nb_aln_done = 0
    progressing = 0

    id_read = ""
    while True:
        line = sam_file.readline()

        if not line:
            break
        if line[0] == "@": # header line, useless
            continue

        # progressing bar
        nb_aln_done += 1
        tmp_progressing = int(nb_aln_done / nb_tot_alignments * 100)
        if tmp_progressing > progressing and tmp_progressing % 1 == 0:
            progressing = tmp_progressing
            display_progressing_bar(nb_aln_done, progressing, time.time())

        line_list = line.split()

        # check if the read is unmapped
        flag = line_list[1]
        if flag in ["4", "256", "2048"]:
            continue
        # check if this is not a secondary alignment
        read = line_list[9]
        if read == "0" or read == "*" or id_read == line_list[0]:
            continue

        id_read = line_list[0]
        genome_position = int(line_list[3])
        cigar = line.split()[5]
        # call function that will get error information and update results:
        parse_cigar(cigar, genome_position, dict_tmp)

    return(dict_tmp)

def display_progressing_bar(nb_aln_done: int, prct: int, current_time: float):
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
    elapsed_time = current_time - time_start
    estimated_remaining_time = (alignment_total_nb-nb_aln_done) * elapsed_time / nb_aln_done
    estimated_remaining_time = str(datetime.timedelta(seconds=round(estimated_remaining_time)))

    # Write to stdout
    sys.stdout.write(f"\r{progressing_bar} ({prct}% ; {estimated_remaining_time})")
    sys.stdout.flush()


def parse_cigar(cigar: str, genome_position: int, dict_tmp: dict):
    """ Takes the CIGAR field as input (string), and turns it into list of elements
        Also updates the dictionary (dict_tmp) storing error rates for each
          relative position in the genome, for the current species
    Parameters
    ----------
    cigar: str
        The parameter `cigar` contains the CIGAR field extracted from SAM file.
        It contains information about number of matches/mismatches, indels
          and soft / hard clips (strings of integers).
    """

    if "reverse" in sam_filename:
        genome_position = reference_genome_length - genome_position
        step_position = -1
    else:
        step_position = 1

    # Parse CIGAR string and split it into a list of elements
    cigar_list = []
    pattern = r"([0-9]*[S=XDIH]{1})"
    for match in re.finditer(pattern, cigar):
        cigar_list += [match.group()]

    for relative_position in range(0, 100):
        if relative_position not in dict_tmp:
            dict_tmp[relative_position] = [0, 0, 0, 0]


    # Updates dict_tmp
    position = genome_position
    for cigar_elt in cigar_list:
        nb_occ = int(cigar_elt[:-1])
        error_type = cigar_elt[-1]
        if error_type in ("S", "H"): # skip: soft/hard clips that are not part of alignment
            continue
        while nb_occ > 0:
            relative_position = int(position / reference_genome_length * 100)
            if relative_position == 100:
                relative_position = 99 # margin effect, very last base of genome
#            if relative_position not in dict_tmp:
#                dict_tmp[relative_position] = [0, 0, 0, 0]
            if error_type == "=": # match
                dict_tmp[relative_position][0] += 1
                position += step_position
            elif error_type == "X": # mismatch
                dict_tmp[relative_position][1] += 1
                position += step_position
            elif error_type == "I": # insertion
                dict_tmp[relative_position][2] += 1
            elif error_type == "D": # deletion
                dict_tmp[relative_position][3] += 1
                position += step_position

            nb_occ -= 1

    return dict_tmp


def sum_lists(list1, list2, nb_occ):
    """
    Sum values of list1 and list2
    """
    if len(list1) != len(list2):
        print("ERROR: lists are not of same length")
        exit(1)
    new_list = [] #[list1[i]+list2[i] for i in range(len(list1))]
    for i in range(len(list1)):
        elt1, elt2 = list1[i], list2[i]
        if math.isnan(elt1) and math.isnan(elt2):
            new_list += [np.nan]
        else:
            new_value = 0
            if not math.isnan(elt1):
                new_value += elt1
                nb_occ[i] += 1
            if not math.isnan(elt2):
                new_value += elt2
                nb_occ[i] += 1
            new_list += [new_value]
        
    return new_list, nb_occ




def compute_results():
    """
    Computes error rates for each relative position, for the current species,
        stored in `dict_error_tmp`,
        and adds results in `dict_error_rates_position`
    Plots the results contained in `dict_error_rates_position` (saved in OUTPUT_PLOT directory)
    """

    # Updates dict_error_rates_position with current species' results
    for position in dict_error_tmp:
        nb_match, nb_mismatch, nb_insertion, nb_deletion = dict_error_tmp[position]
        alignment_total_length = sum([nb_match, nb_mismatch, nb_insertion, nb_deletion])
        if alignment_total_length == 0:
            dict_error_rates_position["Mismatch"][position] += [np.nan]
            dict_error_rates_position["Insertion"][position] += [np.nan]
            dict_error_rates_position["Deletion"][position] += [np.nan]
            dict_error_rates_position["Total"][position] += [np.nan]
            continue
        mismatch_rate = nb_mismatch / alignment_total_length * 100
        insertion_rate = nb_insertion / alignment_total_length * 100
        deletion_rate = nb_deletion / alignment_total_length * 100
        total_error_rate = 100 - (nb_match / alignment_total_length * 100)
        dict_error_rates_position["Mismatch"][position] += [mismatch_rate]
        dict_error_rates_position["Insertion"][position] += [insertion_rate]
        dict_error_rates_position["Deletion"][position] += [deletion_rate]
        dict_error_rates_position["Total"][position] += [total_error_rate]
    dict_mismatches, dict_insertions, dict_deletions, dict_total = {}, {}, {}, {}

    list_error_names = ["Mismatch", "Insertion", "Deletion", "Total"]
    list_dictionaries = [dict_mismatches, dict_insertions, dict_deletions, dict_total]
    for i in range(len(list_error_names)):
        error = list_error_names[i]
        dictionary = list_dictionaries[i]
        for position in sorted(dict_error_rates_position[error]):
            error_rate_list = dict_error_rates_position[error][position]
            for j, error_rate in enumerate(error_rate_list):
                if j >= len(L_species):
                    continue
                species_name = L_species[j]
                if species_name not in dictionary:
                    dictionary[species_name] = []
                dictionary[species_name] += [error_rate]
        list_dictionaries[i] = dictionary
    dict_mismatches, dict_insertions, dict_deletions, dict_total = list_dictionaries
    
    
    # TEST ------

    
    # Plot results - if species are grouped
    if os.path.exists(FILE_SPECIES_GROUPS):
        
        # Get group names and colors
        dict_category_color = {}
        dict_group_species = {}
        color_file = open(FILE_SPECIES_GROUPS, "r")
        for line in color_file:
            group_name, species_name, color = line.rstrip().split("\t")
            species_name = species_name.replace(" ", "_")
            dict_group_species[species_name] = group_name
            dict_category_color[group_name] = color
        color_file.close()
        
        
        # Re-arrange dictionnary before plotting results
        for i in range(len(list_error_names)):
            dictionary = list_dictionaries[i]
            grouped_dictionary = {}
            dict_nb_occ = {}
            for species_name in L_species:
                group_name = dict_group_species[species_name]
                if group_name not in grouped_dictionary:
                    grouped_dictionary[group_name] = dictionary[species_name]
                    nb_occ = [1 if not math.isnan(elt) else 0 for elt in grouped_dictionary[group_name]]
                    dict_nb_occ[group_name] = nb_occ
                else:
                    old_value = grouped_dictionary[group_name]
                    to_add = dictionary[species_name]
                    new_value, nb_occ = sum_lists(old_value, to_add, dict_nb_occ[group_name])
                    dict_nb_occ[group_name] = nb_occ
                    grouped_dictionary[group_name] = new_value
                    
            for group_name in dict_nb_occ:
                nb_occ = dict_nb_occ[group_name]
                value_sum = grouped_dictionary[group_name]
                value_mean = [value_sum[i]/nb_occ[i] for i in range(len(value_sum))]
                grouped_dictionary[group_name] = value_mean
            list_dictionaries[i] = grouped_dictionary
        
        # Plot

        for i in range(len(list_error_names)):
            error = list_error_names[i]
            dictionary = list_dictionaries[i]
    
            fig = plt.figure()
            ax = fig.add_subplot(111)
            fig.set_dpi(300.0)
            box = ax.get_position()
            ax.set_position([0.1, 0.12, box.width*0.8, box.height])
            for group_name in dictionary:
                color = dict_category_color[group_name]
                ax.plot(dictionary[group_name], color=color, label=group_name)
            ax.set(xlabel='Relative position in reference genome (%)',
                   ylabel=error + ' error rate (%)')
    
            # Reorder legend labels
            ax.legend(title="Species groups",
                      fontsize=6, title_fontsize=8,
                       bbox_to_anchor=(1, 1)) # place legend outside plot
            plt.savefig(OUTPUT_PLOT + f"error_rates_along_genome_{error}.png")
            plt.close()
    
    # Plot results - if species are not grouped
    else:
        dict_category_color = {}
        dict_gc_species = {}
        color_file = open(FILENAME_SPECIES_GC_COLOR, "r")
        for line in color_file:
            species_name, color, gc = line.rstrip().split("\t")
            species_name = species_name.replace(" ", "_")
            gc = float(gc)
            dict_category_color[species_name] = color
            dict_gc_species[gc] = species_name
        color_file.close()
        list_ordered_species_gc = []
        for gc in sorted(dict_gc_species):
            list_ordered_species_gc += [dict_gc_species[gc]]

        for i in range(len(list_error_names)):
            error = list_error_names[i]
            dictionary = list_dictionaries[i]
    
            fig = plt.figure()
            ax = fig.add_subplot(111)
            fig.set_dpi(300.0)
            box = ax.get_position()
            ax.set_position([0.1, 0.12, box.width*0.8, box.height])
            for species_name in L_species:
                color = dict_category_color[species_name]
                ax.plot(dictionary[species_name], color=color, label=species_name)
            ax.set(xlabel='Relative position in reference genome (%)',
                   ylabel=error + ' error rate (%)')
    
            # Reorder legend labels
            handles, labels = ax.get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            ordered_label_list = [elt for elt in list_ordered_species_gc if elt in by_label]
            ordered_label_values = [by_label[k] for k in ordered_label_list]
            ordered_label_list_short_name = [get_short_name(elt) for elt in ordered_label_list]
            ax.legend(ordered_label_values, ordered_label_list_short_name, title="Species",
                      fontsize=6, title_fontsize=8,
                       bbox_to_anchor=(1, 1)) # place legend outside plot
            plt.savefig(OUTPUT_PLOT + f"error_rates_along_genome_{error}.png")
            plt.close()

def get_short_name(long_name):
    L_name = long_name.replace("_", " ").split(" ")
    if len(L_name) == 1:
        return long name
    L_name[0] = L_name[0][0] + "."
    short_name = " ".join(L_name)
    return short_name



# --------------------------------------------------------------------------------------------------
# Main

if __name__ == "__main__":


    ## Parses arguments
    NUMBER_EXPECTED_ARGUMENTS = 5
    if len(sys.argv) != NUMBER_EXPECTED_ARGUMENTS + 1:
        print(f"ERROR: Wrong number of arguments: {NUMBER_EXPECTED_ARGUMENTS} expected but {len(sys.argv)-1} given.")
        sys.exit(2)
    SAM_DIRNAME, REF_GEN_DIRNAME, OUTPUT_PLOT, FILENAME_SPECIES_GC_COLOR, FILE_SPECIES_GROUPS = sys.argv[1:]


    # --- Initialize variables for the plot ---
    dict_error_rates_position, L_species = initialize(L_errors_types)


    # --- Processing alignments ---
    for sam_filename in os.listdir(SAM_DIRNAME):

        dict_error_tmp = {} # temporary version of dict_error_rates_position

        # - Processing of forward alignments -
        sam_filepath = SAM_DIRNAME + sam_filename
        is_forward, species_name, reference_genome_filename = checks_is_forward_alignment(sam_filepath)
        if not is_forward:
            continue
        L_species += [species_name]
        alignment_total_nb = get_total_number_of_lines(sam_filepath) # overestimated, because of heading lines
        time_start = time.time()

        sam_file = open(sam_filepath, "r")

        # get reference genome total length
        reference_genome_length = get_reference_genome_length(reference_genome_filename)

        # extract errors from sam mapping file, and update dict_error_tmp
        extract_errors_from_sam(sam_file, alignment_total_nb, dict_error_tmp)

        sam_file.close()
        sys.stdout.write("\n")


        # - Processing of reverse alignments -
        sam_filename = sam_filename.replace(".sam", "_reverse.sam")
        sam_filepath = SAM_DIRNAME + sam_filename
        alignment_total_nb = get_total_number_of_lines(sam_filepath) # overestimated, because of heading lines
        time_start = time.time()

        sam_file = open(sam_filepath, "r")

        # extract errors from sam mapping file, and update dict_error_tmp
        dict_error_tmp = extract_errors_from_sam(sam_file, alignment_total_nb, dict_error_tmp)

        sam_file.close()
        sys.stdout.write("\n")


        # - Compute results -
        compute_results()
        
        
    print("----")
    print(L_species)

    compute_results()




