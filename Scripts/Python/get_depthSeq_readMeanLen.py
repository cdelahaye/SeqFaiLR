#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Computes the sequencing depth, GC content of reference genome, and mean read length
  for all sequenced datasets.
Output results as a table, both in raw file, and in a LaTeX style.
"""

# --------------------------------------------------------------------------------------------------
# Packages

import sys
import os

# --------------------------------------------------------------------------------------------------
# Functions

def get_reference_genome_length_gc(filename) -> int:
    """ Computes the total length of reference genome given as argument
    Expect to be a .fasta file
    """
    genome_length = 0
    nb_g_c_bases = 0
    with open(filename, "r") as ref_genome:
        for line in ref_genome:
            if line[0] != ">": # skip
                line = line.rstrip().upper().replace("N", "")
                genome_length += len(line)
                nb_g_c_bases += line.count("C") + line.count("G")
    gc_content = round(nb_g_c_bases / genome_length * 100, 2)
    return genome_length, gc_content

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
    print("Running.")

    ## Parses arguments
    if len(sys.argv) != 4:
        print(f"ERROR: Wrong number of arguments: 3 expected but {len(sys.argv)-1} given.")
        sys.exit(2)
    ALN_DIRNAME = sys.argv[1]
    REF_GEN_DIRNAME = sys.argv[2]
    OUTPUT_RAW = sys.argv[3]

    if ALN_DIRNAME[-1] != "/":
        ALN_DIRNAME += "/"
    if REF_GEN_DIRNAME[-1] != "/":
        REF_GEN_DIRNAME += "/"
    if OUTPUT_RAW[-1] != "/":
        OUTPUT_RAW += "/"

    dict_species_stats = {}

    # Browse each alignment file
    for aln_filename in os.listdir(ALN_DIRNAME):
        if not os.path.isfile(ALN_DIRNAME + aln_filename):
            continue

        print(aln_filename)

        species_name = aln_filename.split(".")[0]

        # Compute reference genome length
        reference_filename = REF_GEN_DIRNAME + species_name + ".fasta"
        ref_genome_length, ref_genome_gc = get_reference_genome_length_gc(reference_filename)

        # Parse alignment file to get coverage and mean read length for aligned reads
        sum_aligned_read_bases = 0
        nb_aligned_read = 0

        aln_file = open(ALN_DIRNAME + aln_filename, "r")

        while True:
            aln_file.readline() # header, skip
            aln_file.readline() # genome, skip
            read = aln_file.readline().replace("-", "")
            if not read:
                break
            read_len = len(read)
            sum_aligned_read_bases += read_len
            nb_aligned_read += 1

        aln_file.close()

        depth_seq = str(round(sum_aligned_read_bases / ref_genome_length)) + "X"
        mean_read_len = int(round(sum_aligned_read_bases / nb_aligned_read))
        species_name_short = get_short_name(species_name)
        dict_species_stats[species_name_short] = [ref_genome_length, ref_genome_gc,
                                                  mean_read_len, depth_seq]



    # Store results in file

    file = open(OUTPUT_RAW + "sequencing_depth_mean_read_length.txt", "w")

    # Tabulated format
    header_list = ["Species", "Size (non-N bases)", "GC %",
                   "Mean read length (bp)", "Sequencing depth"]
    file.write("\t".join(header_list) + "\n")
    for species_name in sorted(dict_species_stats):
        file.write(species_name + "\t")
        file.write("\t".join(list(map(str, dict_species_stats[species_name]))) + "\n")

    # LaTeX format
    file.write("\n")
    file.write("\\begin{table}[!ht]\n\centering\n\\begin{tabular}{lcccc}\n")
    header_list = ["Species", "Size (non-N bases)", "GC \%",
                   "Mean read length (bp)", "Sequencing depth"]
    file.write(" & ".join(header_list) + " \\\\ \n")
    for species_name in sorted(dict_species_stats):
        file.write("\\textit{" + species_name + "} & ")
        file.write(" & ".join(list(map(str, dict_species_stats[species_name]))) + " \\\\ \n")
    file.write("\end{tabular} \n\end{table}")


    file.close()
