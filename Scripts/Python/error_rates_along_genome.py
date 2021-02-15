#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Computes error rates (mismatches, insertions, deletions, and total) depending on the relative
position in the reference genome.
"""

# --------------------------------------------------------------------------------------------------
# Packages

import re
import sys
import os
import time
import datetime
import numpy as np
import matplotlib.pyplot as plt
# Parameters for output plots
PARAMS = {'legend.fontsize': 18,
          'legend.title_fontsize': 18,
          'legend.labelspacing':0.1,
          'legend.borderpad':0.3,
          'legend.columnspacing':1,
          'legend.handletextpad':0.5,
          'legend.handlelength':0.8,
          'figure.figsize': (14, 9),
          'axes.labelsize': 20,
          'axes.titlesize':22,
          'xtick.labelsize':20,
          'ytick.labelsize':20}
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


def get_reference_genome_length(filename) -> int:
    """ Computes the total length of reference genome given as argument
    Expect to be a .fasta file
    """
    genome_length = 0
    with open(filename, "r") as ref_genome:
        for line in ref_genome:
            if line[0] != ">": # skip
                genome_length += len(line.replace("\n", ""))
    return genome_length


def parse_cigar(cigar: str):
    """ Takes the CIGAR field as input (string), and turns it into list of elements
        Also updates the dictionary (location_error_rates_dict_temp) storing error rates for each
          relative position in the genome, for the current species
    Parameters
    ----------
    cigar: str
        The parameter `cigar` contains the CIGAR field extracted from SAM file.
        It contains information about number of matches/mismatches, indels
          and soft / hard clips (strings of integers).
    """

    if "reverse" in sam_filename:
        step_position = -1
    else:
        step_position = 1

    # Parse CIGAR string and split it into a list of elements
    cigar_list = []
    pattern = r"([0-9]*[S=XDIH]{1})"
    for match in re.finditer(pattern, cigar):
        cigar_list += [match.group()]

    # Updates location_error_rates_dict_temp
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
            if relative_position not in location_error_rates_dict_temp:
                location_error_rates_dict_temp[relative_position] = [0, 0, 0, 0]
            if error_type == "=": # match
                location_error_rates_dict_temp[relative_position][0] += 1
                position += step_position
            elif error_type == "X": # mismatch
                location_error_rates_dict_temp[relative_position][1] += 1
                position += step_position
            elif error_type == "I": # insertion
                location_error_rates_dict_temp[relative_position][2] += 1
            elif error_type == "D": # deletion
                location_error_rates_dict_temp[relative_position][3] += 1
                position += step_position

            nb_occ -= 1


def compute_results():
    """Computes error rates for each relative position, for the current species,
        stored in `location_error_rates_dict_temp`,
        and adds results in `location_error_rates_dict_toplot`
    Plots the results contained in `location_error_rates_dict_toplot` (saved in OUTPUT_PLOT directory),
    and write raw results (in OUTPUT_RAW directory)
    """

    # Updates location_error_rates_dict_toplot with current species' results
    for position in location_error_rates_dict_temp:
        nb_match, nb_mismatch, nb_insertion, nb_deletion = location_error_rates_dict_temp[position]
        alignment_total_length = sum([nb_match, nb_mismatch, nb_insertion, nb_deletion])
        mismatch_rate = nb_mismatch / alignment_total_length * 100
        insertion_rate = nb_insertion / alignment_total_length * 100
        deletion_rate = nb_deletion / alignment_total_length * 100
        total_error_rate = 100 - (nb_match / alignment_total_length * 100)
        location_error_rates_dict_toplot["Mismatch"][position] += [mismatch_rate]
        location_error_rates_dict_toplot["Insertion"][position] += [insertion_rate]
        location_error_rates_dict_toplot["Deletion"][position] += [deletion_rate]
        location_error_rates_dict_toplot["Total"][position] += [total_error_rate]

    # Plot lines: error rate for each species
    dict_mismatches, dict_insertions, dict_deletions, dict_total = {}, {}, {}, {}
    for position in sorted(location_error_rates_dict_toplot["Mismatch"]):
        mismatch_rate_list = location_error_rates_dict_toplot["Mismatch"][position]
        for i, mismatch_rate in enumerate(mismatch_rate_list):
            if i not in dict_mismatches:
                dict_mismatches[i] = []
            dict_mismatches[i] += [mismatch_rate]
    for position in sorted(location_error_rates_dict_toplot["Insertion"]):
        insertion_rate_list = location_error_rates_dict_toplot["Insertion"][position]
        for i, insertion_rate in enumerate(insertion_rate_list):
            if i not in dict_insertions:
                dict_insertions[i] = []
            dict_insertions[i] += [insertion_rate]
    for position in sorted(location_error_rates_dict_toplot["Deletion"]):
        deletion_rate_list = location_error_rates_dict_toplot["Deletion"][position]
        for i, deletion_rate in enumerate(deletion_rate_list):
            if i not in dict_deletions:
                dict_deletions[i] = []
            dict_deletions[i] += [deletion_rate]
    for position in sorted(location_error_rates_dict_toplot["Total"]):
        total_error_rate_list = location_error_rates_dict_toplot["Total"][position]
        for i, total_error_rate in enumerate(total_error_rate_list):
            if i not in dict_total:
                dict_total[i] = []
            dict_total[i] += [total_error_rate]
    fig, axs = plt.subplots(2, 2)
    plt.subplots_adjust(left=0.05, right=0.95, bottom=0.1, top=0.95,
                        wspace=0.2, hspace=0.4)
    fig.set_dpi(300.0)
    for species_id in range(len(species_list)):
        axs[0, 0].plot(dict_mismatches[species_id], color=COLOR_MISMATCH)
        axs[0, 0].set_title('A. Mismatches')
        axs[0, 1].plot(dict_insertions[species_id], color=COLOR_INSERTION)
        axs[0, 1].set_title('B. Insertions')
        axs[1, 0].plot(dict_deletions[species_id], color=COLOR_DELETION)
        axs[1, 0].set_title('C. Deletions')
        axs[1, 1].plot(dict_total[species_id], color="k")
        axs[1, 1].set_title('D. Total')
    for ax in axs.flat:
        ax.set(xlabel='Relative position in reference genome (%)',
               ylabel='Error rate (%)',
               ylim=(0, 6))
        ax.set_yticks(np.arange(0, 7))
        ax.set_yticklabels(np.arange(0, 7))
    axs[1, 1].set(ylim=(2, 10))
    axs[1, 1].set_yticks(np.arange(2, 11))
    axs[1, 1].set_yticklabels(np.arange(2, 11))
    axs[0, 0].set(xlabel="")
    axs[0, 1].set(xlabel="", ylabel="")
    axs[1, 1].set(ylabel="")
    plt.savefig(OUTPUT_PLOT + "error_rates_along_genome_lines.png")
    plt.close()

    # Same plot, only keep one of the four plots: deletions

    dict_species_color = {}
    dict_gc_species = {}
    color_file = open(FILENAME_SPECIES_COLOR, "r")
    for line in color_file:
        species_name, color, gc = line.rstrip().split(" ; ")
        species_name = species_name.replace(" ", "_")
        gc = float(gc)
        if color[0] != "#":
            color = "#" + color
        dict_species_color[species_name] = color
        dict_gc_species[gc] = species_name
    color_file.close()
    list_ordered_species_gc = []
    for gc in sorted(dict_gc_species):
        list_ordered_species_gc += [get_short_name(dict_gc_species[gc])]
    # plot
    fig, ax = plt.subplots()
    fig.set_dpi(300.0)
    for i, species_name in enumerate(species_list):
        color = dict_species_color[species_name]
        ax.plot(dict_deletions[i], color=color, label=get_short_name(species_name))
    ax.set(xlabel='Relative position in reference genome (%)',
           ylabel='Deletions error rate (%)',
           ylim=(0, 6))
    ax.set_yticks(np.arange(0, 6, 1))
    ax.set_yticklabels(np.arange(0, 6, 1))

    # Reorder legend labels
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ordered_label_list = [elt for elt in list_ordered_species_gc if elt in by_label]
    ordered_label_values = [by_label[k] for k in ordered_label_list]
    plt.legend(ordered_label_values, ordered_label_list, ncol=3, title="Species")


    plt.savefig(OUTPUT_PLOT + "error_rates_along_genome_lines_deletions_only.png")
    plt.close()


def get_short_name(long_species_name: str):
    """ Returns a short version of species name
         and replace underscores with spaces
    For example: Staphylococcus_thermophilus_CNRZ1066 -> S. thermophilus CNRZ1066
    """
    list_long_species_name = long_species_name.split("_")
    short_species_name = list_long_species_name[0][0] + ". "
    short_species_name += " ".join(list_long_species_name[1:])
    return short_species_name


# --------------------------------------------------------------------------------------------------
# Main

if __name__ == "__main__":

    # Colors for plot:
    COLOR_MISMATCH = "#DDAA33"
    COLOR_INSERTION = "#BB5566"
    COLOR_DELETION = "#004488"
    L_COLORS = [COLOR_MISMATCH, COLOR_INSERTION, COLOR_DELETION, "k"]

    ## Parses arguments

    if len(sys.argv) != 6:
        print(f"ERROR: Wrong number of arguments: 5 expected but {len(sys.argv)-1} given.")
        sys.exit(2)
    SAM_DIRNAME = sys.argv[1]
    REF_GEN_DIRNAME = sys.argv[2]
    OUTPUT_RAW = sys.argv[3]
    OUTPUT_PLOT = sys.argv[4]
    FILENAME_SPECIES_COLOR = sys.argv[5]

    if SAM_DIRNAME[-1] != "/":
        SAM_DIRNAME += "/"
    if REF_GEN_DIRNAME[-1] != "/":
        REF_GEN_DIRNAME += "/"
    if OUTPUT_PLOT[-1] != "/":
        OUTPUT_PLOT += "/"
    if OUTPUT_RAW[-1] != "/":
        OUTPUT_RAW += "/"

    location_error_rates_dict_toplot = {}
    for error_type in ["Mismatch", "Insertion", "Deletion", "Total"]:
        location_error_rates_dict_toplot[error_type] = {}
        for position in range(100):
            location_error_rates_dict_toplot[error_type][position] = []
    species_list = []

    for sam_filename in os.listdir(SAM_DIRNAME):

        if not os.path.isfile(SAM_DIRNAME + sam_filename):
            continue

        # Reverse mapping files are treated right after their forward homologous file
        # they are not ignored, but not handled here
        if "reverse" in sam_filename:
            continue

        species_name = sam_filename.split(".sam")[0]

        # Ignore datasetss for which genomes are in multiple chromosomes ; only those with a single one are handled
        REF_GEN_FILENAME = REF_GEN_DIRNAME + species_name + ".fasta"
        count_start_fasta = 0
        file = open(REF_GEN_FILENAME, "r")
        for line in file:
            if line[0] == ">":
                count_start_fasta += 1
            if count_start_fasta > 1:
                break
        file.close()
        if count_start_fasta > 1:
            print(f"{species_name} ignored because its reference genome is in multiple part")
            continue


        # stores error rates for relative position in genome, for the current species
        #  of form: dictionary[genomic_position] = [#match, #mismatch, #insertion, #deletion]
        location_error_rates_dict_temp = {}


        ## 1/ Processing of forward alignments
        sam_path = SAM_DIRNAME + sam_filename
        NB_TOT_ALN = get_total_number_of_lines(sam_path) # overestimated, because of heading lines
        STARTING_TIME = time.time()

        print(sam_filename)
        sam_file = open(sam_path, "r")
        species_list += [species_name]

        # Get reference genome total length
        print("Computing reference genome length.")
        reference_genome_length = get_reference_genome_length(REF_GEN_FILENAME)
        print("Done.")


        # Exctract errors from sam mapping file
        nb_aln_done = 0
        progressing = 0
        id_read = ""
        while True:
            line = sam_file.readline()
            if not line:
                break
            if line[0] == "@": # header line, useless
                continue
            nb_aln_done += 1
            tmp_progressing = int(nb_aln_done / NB_TOT_ALN * 100)
            if tmp_progressing > progressing and tmp_progressing % 1 == 0:
                progressing = tmp_progressing
                display_progressing_bar(progressing, time.time())

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
            parse_cigar(cigar)

        sam_file.close()
        sys.stdout.write("\n")


        ## 2/ Processing of reverse alignments
        sam_filename = sam_filename.replace(".sam", "_reverse.sam")
        sam_path = SAM_DIRNAME + sam_filename
        NB_TOT_ALN = get_total_number_of_lines(sam_path) # overestimated, because of heading lines
        STARTING_TIME = time.time()

        print(sam_filename)
        sam_file = open(sam_path, "r")

        # Exctract errors from sam mapping file
        nb_aln_done = 0
        progressing = 0
        id_read = ""
        while True:
            line = sam_file.readline()
            if not line:
                break
            if line[0] == "@": # header line, useless
                continue
            nb_aln_done += 1
            tmp_progressing = int(nb_aln_done / NB_TOT_ALN * 100)
            if tmp_progressing > progressing and tmp_progressing % 1 == 0:
                progressing = tmp_progressing
                display_progressing_bar(progressing, time.time())

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
            genome_position = reference_genome_length - int(line_list[3])
            cigar = line.split()[5]
            parse_cigar(cigar)

        sam_file.close()
        sys.stdout.write("\n")


        compute_results()

    compute_results()
