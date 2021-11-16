#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
 As alignments of reads on forward and reverse strand has been computed separately,
  a given read may have been aligned on both
 Here we parse sam file and keep the best alignment for each read
 Based on soft clip length, because a read mapping on a given strand will have poor
  alignment on the other one, i.e. high portion of soft clipped read
 + Also removing alignments for which the soft clipped length is greater than one half
    of the initial read length
"""

# --------------------------------------------------------------------------------------------------
# Packages

import sys
import ast
import re

# --------------------------------------------------------------------------------------------------
# Functions

def get_soft_clips(cigar: str):
    """ Parse the cigar and returns the soft clips

    Parameters
    ----------
    cigar: str
        The parameter `cigar` contains the CIGAR field extracted from SAM file.
        It contains information about number of matches/mismatches, indels
          and soft / hard clips (strings of integers).

    Returns
    -------
    soft_clip_start, soft_clip_end: int
        number of bases soft-clipped
    """

    # Parse CIGAR string and split it into a list of elements
    cigar_list = []
    pattern = r"([0-9]*[S=XDIH]{1})"
    for match in re.finditer(pattern, cigar):
        cigar_list += [match.group()]

    # Hard clips represent non aligned parts of the reads (3' and 5')
    #  that have already been removed from read sequence
    #  so we just remove hard clip elements from CIGAR field
    if "H" in cigar_list[0]:
        cigar_list = cigar_list[1:]
    if "H" in cigar_list[-1]:
        cigar_list = cigar_list[:-1]

    # Soft clips represent non aligned parts of the reads (3' and 5')
    #  that have to be removed from read sequence
    nb_clipped_start, nb_clipped_end = 0, 0
    if "S" in cigar_list[0]:
        soft_clip = cigar_list[0]
        nb_clipped_start = int(soft_clip[:-1])
    if "S" in cigar_list[-1]:
        soft_clip = cigar_list[-1]
        nb_clipped_end = int(soft_clip[:-1])

    return([nb_clipped_start, nb_clipped_end])


def get_information(SAM_FILENAME):
    """ Parses the given file and store soft clips for each read
    """
    read_id = ""
    file = open(SAM_FILENAME, "r")
    for line in file:
        if line[0] == "@":
            continue
        line = line.rstrip()
        list_line = line.split("\t")

        # check if the read is unmapped
        flag = list_line[1]
        if flag in ["4", "256", "2048"]:
            continue
        # check if this is not a secondary alignment
        read = list_line[9]
        if read == "0" or read == "*" or read_id == list_line[0]:
            continue

        read_id = list_line[0]
        cigar_str = list_line[5]
        list_soft_clips = get_soft_clips(cigar_str)

        if read_id not in dict_read_softclips:
            dict_read_softclips[read_id] = {}
        if strand in dict_read_softclips[read_id]:
            raise Exception(f"The read {read_id} has been mapped (at least) twice on strand",
                            f"{strand}, this should not happen")
        dict_read_softclips[read_id][strand] = list_soft_clips
    file.close()


def delete_wrong_alignments(SAM_FILENAME):
    """ Re-write sam file without wrong alignments (as defined in begining of this python file)
    """
    file = open(SAM_FILENAME, "r")
    new_file = open(SAM_FILENAME.replace(".sam", "_clean.sam"), "w")
    other_strand = {"+": "-", "-": "+"}[strand]
    read_id = ""
    for line in file:
        if line[0] == "@":
            continue
        line = line.rstrip()
        list_line = line.split("\t")

        # check if the read is unmapped
        flag = list_line[1]
        if flag in ["4", "256", "2048"]:
            continue
        # check if this is not a secondary alignment
        read = list_line[9]
        if read == "0" or read == "*" or read_id == list_line[0]:
            continue

        read_id = list_line[0]
        read_length = len(list_line[9])

        to_keep = False # boolean to know if a read should be kept or not

        if len(dict_read_softclips[read_id]) == 1 and strand in dict_read_softclips[read_id]:
            to_keep = True
        elif len(dict_read_softclips[read_id]) == 2:
            sum_soft_clip_strand = sum(dict_read_softclips[read_id][strand])
            sum_soft_clip_other_strand = sum(dict_read_softclips[read_id][other_strand])
            diff_soft_clip = sum_soft_clip_strand - sum_soft_clip_other_strand

            if diff_soft_clip < 0:
                to_keep = True
            if sum_soft_clip_strand / read_length > 0.5:
                to_keep = False

        if not to_keep:
            continue

        new_file.write(line + "\n")

    file.close()
    new_file.close()




# --------------------------------------------------------------------------------------------------
# Main

if __name__ == "__main__":


    ## Parses arguments
    if len(sys.argv) != 3:
        print(f"ERROR: Wrong number of arguments: 2 expected but {len(sys.argv)-1} given.")
        sys.exit(2)
    SPECIES_NAME = sys.argv[1]
    SAM_DIR = sys.argv[2]

    if SAM_DIR[-1] != "/":
        SAM_DIR += "/"

    # Dictionary that will store soft clipped length for each read and each strand
    dict_read_softclips = {}

    # Retrieve information for forward strand
    SAM_FILE_FORWARD = SAM_DIR + SPECIES_NAME + ".sam"
    strand = "+"
    get_information(SAM_FILE_FORWARD)

    # Retrieve information for reverse strand
    SAM_FILE_REVERSE = SAM_DIR + SPECIES_NAME + "_reverse.sam"
    strand = "-"
    get_information(SAM_FILE_REVERSE)


    # Delete bad forward alignments
    strand = "+"
    delete_wrong_alignments(SAM_FILE_FORWARD)


    # Delete bad forward alignments
    strand = "-"
    delete_wrong_alignments(SAM_FILE_REVERSE)

