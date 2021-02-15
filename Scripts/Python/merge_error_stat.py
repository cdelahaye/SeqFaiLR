#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
When explicit alignment (genome versus read) is computed from .sam file, a file is generated
  containing some figures about sequencing errors
As explicit alignment for each species is computed separately for forward and reverse mapping,
  this script aims at merging the two output files into a single one
"""

# --------------------------------------------------------------------------------------------------
# Packages

import sys

# --------------------------------------------------------------------------------------------------
# Functions

def parse_file(filename, total_aln_length, nb_mismatch, nb_insertion, nb_deletion, nb_total_errors):
    """ Parse a given file to extract its information ; updates global dictionaries and variables
    """
    file = open(filename, "r")

    total_aln_length += int(file.readline().rstrip().split(" ")[-1])
    file.readline() # skip global error rate %
    nb_insertion += float(file.readline().rstrip().split(" ")[-1])
    nb_deletion += float(file.readline().rstrip().split(" ")[-1])

    file.readline() # skip empty line
    file.readline() # skip header
    file.readline() # skip another header
    line = file.readline().rstrip()
    while line != "":
        substi, nb_occurrences = line.split("\t")
        nb_occurrences = int(nb_occurrences)
        if substi not in dict_substitution:
            dict_substitution[substi] = 0
        dict_substitution[substi] += nb_occurrences
        nb_mismatch += nb_occurrences
        nb_total_errors += nb_occurrences
        line = file.readline().rstrip()

    file.readline() # skip header
    file.readline() # skip another header
    line = file.readline().rstrip()
    while line != "":
        insertion_length, nb_occurrences = line.split("\t")
        insertion_length = int(insertion_length)
        nb_occurrences = int(nb_occurrences)
        if insertion_length not in dict_insertion_length:
            dict_insertion_length[insertion_length] = 0
        dict_insertion_length[insertion_length] += nb_occurrences
        nb_insertion += nb_occurrences * insertion_length
        nb_total_errors += nb_occurrences * insertion_length
        line = file.readline().rstrip()

    file.readline() # skip header
    file.readline() # skip another header
    line = file.readline().rstrip()
    while line != "":
        deletion_length, nb_occurrences = line.split("\t")
        deletion_length = int(deletion_length)
        nb_occurrences = int(nb_occurrences)
        if deletion_length not in dict_deletion_length:
            dict_deletion_length[deletion_length] = 0
        dict_deletion_length[deletion_length] += nb_occurrences
        nb_deletion += nb_occurrences * deletion_length
        nb_total_errors += nb_occurrences * deletion_length
        line = file.readline().rstrip()

    file.close()

    return total_aln_length, nb_mismatch, nb_insertion, nb_deletion, nb_total_errors


# --------------------------------------------------------------------------------------------------
# Main

if __name__ == "__main__":

    # Parsing arguments
    if len(sys.argv) != 4:
        print(f"ERROR: Wrong number of arguments: 3 expected but {len(sys.argv)-1} given.")
        sys.exit(2)
    FILENAME_STAT_ERR_FORWARD = sys.argv[1]
    FILENAME_STAT_ERR_REVERSE = sys.argv[2]
    OUTPUT_FILENAME = sys.argv[3]

    # Initiate dictionaries and variables
    total_aln_length = 0
    nb_mismatch, nb_insertion, nb_deletion, nb_total_errors = 0, 0, 0, 0
    dict_substitution = {}
    dict_insertion_length = {}
    dict_deletion_length = {}


    # Parse forward file
    total_aln_length, nb_mismatch, nb_insertion, nb_deletion, nb_total_errors = parse_file(FILENAME_STAT_ERR_FORWARD,
                                                                                           total_aln_length,
                                                                                           nb_mismatch,
                                                                                           nb_insertion,
                                                                                           nb_deletion,
                                                                                           nb_total_errors)
    # Parse reverse file
    total_aln_length, nb_mismatch, nb_insertion, nb_deletion, nb_total_errors = parse_file(FILENAME_STAT_ERR_REVERSE,
                                                                                           total_aln_length,
                                                                                           nb_mismatch,
                                                                                           nb_insertion,
                                                                                           nb_deletion,
                                                                                           nb_total_errors)

    # Write gathered results in output file
    file = open(OUTPUT_FILENAME, "w")

    file.write(f"Total alignment length: {total_aln_length}\n")
    file.write(f"Global error rate (%): {round(nb_total_errors / total_aln_length * 100, 2)}\n")
    file.write(f"Insertion rate (%): {round(nb_insertion / total_aln_length * 100, 2)}\n")
    file.write(f"Deletion rate (%): {round(nb_deletion / total_aln_length * 100, 2)}\n\n")

    file.write("Substitutions (genome to read):\nSubstitutions\tOccurrences\n")
    for substitution in sorted(dict_substitution):
        if substitution == "N":
            continue
        file.write(f"{substitution}\t{dict_substitution[substitution]}\n")
    file.write(f"N\t{dict_substitution['N']}\n\n")

    file.write("Insertion lengths distributions:\nLength\tOccurrences\n")
    for insertion_length in sorted(dict_insertion_length):
        file.write(f"{insertion_length}\t{dict_insertion_length[insertion_length]}\n")
    file.write("\n")

    file.write("Deletion lengths distributions:\nLength\tOccurrences\n")
    for deletion_length in sorted(dict_deletion_length):
        file.write(f"{deletion_length}\t{dict_deletion_length[deletion_length]}\n")
    file.write("\n")

    file.close()
