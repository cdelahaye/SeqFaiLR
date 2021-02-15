#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program retrieves alignment of reads against a genome, by using MD field and CIGAR string
  provided in the input .sam file
It outputs the final alignments with, for each read, its alignment against the genome
  (the genome is first printed, then the read)
"""

# --------------------------------------------------------------------------------------------------
# Packages

import re
import sys
import time
import datetime


# --------------------------------------------------------------------------------------------------
# Classes

class Alignment:
    """ Explicit alignment between a read and a genome part.
    Attributes:
    -----------
    genome_aln: str
        Alignment of the genome
    read_aln: str
        Alignment of the read
    cigar: [str]
        The parameter `cigar` contains the CIGAR field extracted from SAM file.
        It contains information about number of matches/mismatches, indels
          and soft / hard clips (strings of integers).
    md: [str]
        The parameter `md` contains the MD field extracted from SAM file.
        It contains information about: matches/mismatches combined (integers),
        insertions (letters), and deletions (letters preceded by ^).
    """

    def __init__(self, read: str, cigar: str, md: str):
        self.read = read # unaligned read
        self.genome_aln = ""
        self.read_aln = ""
        self.cigar = cigar
        self.md = md
        self.match = 0 # number of matches in alignment
        self.mismatch = 0
        self.insertion = 0
        self.deletion = 0
        self.insertion_length_dict = {} # dictionary that will store insertion length occurrences
        self.deletion_length_dict = {}
        self.substitution_counts_dict = {} # Store occurrences of each substitution


    def add_match(self, match_length: int):
        """ Add match of given length to alignment
        """
        self.genome_aln += self.read[: match_length]
        self.read_aln += self.read[: match_length]
        self.read = self.read[match_length:]
        self.match += match_length


    def add_mismatch(self, mismatch_genome_part: str):
        """ Add mismatch, with given bases in genome
        """
        mismatch_length = len(mismatch_genome_part)
        mismatch_read_part = self.read[: mismatch_length]
        self.genome_aln += mismatch_genome_part
        self.read_aln += mismatch_read_part
        self.read = self.read[mismatch_length:]
        self.mismatch += mismatch_length
        for i in range(mismatch_length):
            if mismatch_genome_part[i] == "N":
                substitution = "N"
            else:
                substitution = mismatch_genome_part[i] + mismatch_read_part[i]
            if substitution not in self.substitution_counts_dict:
                self.substitution_counts_dict[substitution] = 0
            self.substitution_counts_dict[substitution] += 1


    def add_insertion(self, insertion_length: int):
        """ Add insertion of given length to alignment
        """
        self.genome_aln += "-" * insertion_length
        self.read_aln += self.read[: insertion_length]
        self.read = self.read[insertion_length:]
        self.insertion += insertion_length
        if insertion_length not in self.insertion_length_dict:
            self.insertion_length_dict[insertion_length] = 0
        self.insertion_length_dict[insertion_length] += 1


    def add_deletion(self, deleted_genome_part: str):
        """ Add deletion, with given genome bases, in alignment
        """
        deleted_genome_part = deleted_genome_part[1:] # remove ^
        self.genome_aln += deleted_genome_part
        deletion_length = len(deleted_genome_part)
        self.read_aln += "-" * deletion_length
        self.deletion += deletion_length
        if deletion_length not in self.deletion_length_dict:
            self.deletion_length_dict[deletion_length] = 0
        self.deletion_length_dict[deletion_length] += 1


    def compute_explicit_alignment(self):
        """ Computes explicit alignment for given read, with help of CIGAR and MD information
        Parameters
        ----------
        self: class Alignment
            Instance of Alignment class
        Returns
        -------
        self:
            Computed alignement and some statistics about it (error rates,
            indel length distribution, ...)
        """

        pos_in_cigar = 0
        while pos_in_cigar < len(self.cigar):
            elt_cigar = self.cigar[pos_in_cigar]

            # Match (and insertions)
            if "=" in elt_cigar:
                elt_md = self.md[: len(elt_cigar)-1]
                self.md = self.md[len(elt_cigar)-1:]
                while self.md != "" and self.md[0].isdigit():
                    elt_md += self.md[0]
                    self.md = self.md[1:]

                nb_match_md = int(elt_md)
                nb_match_cigar = int(elt_cigar[:-1])
                self.add_match(nb_match_cigar)

                # Consecutive matches and insertions in CIGAR are sometimes merged as
                #  matches in MD:
                while nb_match_cigar < nb_match_md:
                    pos_in_cigar += 1
                    elt_cigar = self.cigar[pos_in_cigar]

                    if "I" in elt_cigar:
                        insertion_length = int(elt_cigar[:-1])
                        self.add_insertion(insertion_length)
                    elif "=" in elt_cigar:
                        match_length = int(elt_cigar[:-1])
                        nb_match_cigar += match_length
                        self.add_match(match_length)
                    else:
                        raise ValueError(f"Expected to found match or insertion,\
                                          but found {elt_cigar} instead")
                if nb_match_cigar != nb_match_md:
                    raise ValueError(f"Final length of match of CIGAR ({nb_match_cigar}) differs\
                                     from the one in MD ({nb_match_md}).¯")

            # Mismatch
            elif "X" in elt_cigar:
                nb_mismatch = int(elt_cigar[:-1])
                mismatch_genome_part = self.md[: nb_mismatch * 2 - 1]
                self.md = self.md[nb_mismatch * 2 - 1:]
                if self.md[0] == "0":
                    self.md = self.md[1:]
                if mismatch_genome_part[0] == "0":
                    mismatch_genome_part += self.md[0]
                    self.md = self.md[1:]
                mismatch_genome_part = mismatch_genome_part.replace("0", "")
                if len(mismatch_genome_part) != nb_mismatch:
                    raise ValueError(f"Mismatch length in CIGAR ({nb_mismatch}) differs from the\
                                     one in MD ({len(mismatch_genome_part)}).")
                self.add_mismatch(mismatch_genome_part)


            # Insertion
            elif "I" in elt_cigar:
                insertion_length = int(elt_cigar[:-1])
                self.add_insertion(insertion_length)

            # Deletion
            elif "D" in elt_cigar:
                nb_deletion = int(elt_cigar[:-1])
                deleted_genome_part = self.md[: nb_deletion + 1]
                self.md = self.md[nb_deletion + 1:]
                self.add_deletion(deleted_genome_part)

            # Move on to next CIGAR element
            pos_in_cigar += 1



    def checks_consistency(self):
        """ Ensures that the computed alignment is consistent
        """

        len_genome_aln = len(self.genome_aln)
        len_read_aln = len(self.read_aln)
        sum_errors = self.mismatch + self.insertion + self.deletion

        sum_insertions = 0
        for insertion_length in self.insertion_length_dict:
            sum_insertions += insertion_length * self.insertion_length_dict[insertion_length]

        sum_deletions = 0
        for deletions_length in self.deletion_length_dict:
            sum_deletions += deletions_length * self.deletion_length_dict[deletions_length]

        if len_genome_aln != len_read_aln:
            raise Exception(f"Alignment length of genome ({len_genome_aln}) differs\
                            from alignment length of the read ({len_read_aln}).")

        if self.match + sum_errors != len_genome_aln:
            raise Exception(f"Sum of sequencing errors ({sum_errors}) and \
                            matches ({self.match}) differs from total alignment length\
                            ({len_genome_aln}).")

        if self.insertion != sum_insertions:
            raise Exception(f"Sum of insertion errors in dictionary insertion_length_dict \
                            ({sum_insertions}) differs from total insertions counted\
                            ({self.insertion}).")

        if self.deletion != sum_deletions:
            raise Exception(f"Sum of deletion errors in dictionary deletion_length_dict \
                            ({sum_deletions}) differs from total deletions counted\
                            ({self.deletion}).")

# --------------------------------------------------------------------------------------------------
# Functions

def parse_cigar(cigar: str, read: str, pos_in_g: int):
    """ Takes the CIGAR field as input (string), and turns it into list of elements
        Also update the read sequence and its position in genome,
          taking soft/hard clips into account
    Parameters
    ----------
    cigar: str
        The parameter `cigar` contains the CIGAR field extracted from SAM file.
        It contains information about number of matches/mismatches, indels
          and soft / hard clips (strings of integers).
    read: str
        Read sequence that may have to be shorten because of soft/hard clips
    pos_in_g: str
        Initial position in genome where the read have been mapped initially
        This number can be modified due to soft clips
    Returns
    -------
    cigar_list: [str]
        list of elements of cigar
    read and pos_in_g:
        updated values of input data, taking soft and hard clips into account
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
        cigar_list = cigar_list[1:]
        nb_clipped_start = int(soft_clip[:-1])
        read = read[nb_clipped_start:]
    if "S" in cigar_list[-1]:
        soft_clip = cigar_list[-1]
        cigar_list = cigar_list[:-1]
        nb_clipped_end = int(soft_clip[:-1])
        read = read[:-nb_clipped_end]

    return(cigar_list, read, str(pos_in_g), [nb_clipped_start, nb_clipped_end])


def initialize_substitution_dictionary() -> dict:
    """ Initializes a dictionary that will store occurrences for each possible substitution error
    Returns
    -------
    substitution_occurrences_dict: dict
        keys: all possible substitutions within DNA alphabet (+ if "N" is in reference genome)
        values: occurrences of the associated substitution
    """

    alphabet = ["A", "C", "G", "T"]
    substi_occ_dict = {}
    for base_reference in alphabet:
        for base_read in alphabet:
            if base_reference != base_read:
                substitution = base_reference + base_read
                substi_occ_dict[substitution] = 0
    substi_occ_dict["N"] = 0 # in case an N is in reference genome
    return substi_occ_dict


def update_dictionary(dict_to_add: dict, dict_to_update: dict) -> dict:
    """ Update a given dictionary by adding values of another dictionary
    Will only works for dictionaries of form: dict[key] = int
    Parameters
    ----------
    dict_to_add: dict
        Dictionary containing values we want to add to the dictionary to update
    dict_to_update: dict
        Dictionary we want to update
    Returns
    -------
    dict_to_update: dict
        Updated dictionary
    """
    for key in dict_to_add:
        value = dict_to_add[key]
        if not isinstance(value, int): # Checks structure of the dictionary dict_to_add
            raise TypeError(f"Dictionary not of the expected form")
        if key not in dict_to_update:
            dict_to_update[key] = 0
        dict_to_update[key] += value
    return dict_to_update


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
    estimated_remaining_time = (NB_ALN_TO_COMPUTE-nb_aln_computed) * elapsed_time / nb_aln_computed
    estimated_remaining_time = str(datetime.timedelta(seconds=round(estimated_remaining_time)))

    # Write to stdout
    sys.stdout.write(f"\r{progressing_bar} ({prct}% ; {estimated_remaining_time})")
    sys.stdout.flush()



# --------------------------------------------------------------------------------------------------
# Main

if __name__ == "__main__":

    ## Parses arguments
    if len(sys.argv) != 5:
        print(f"ERROR: Wrong number of arguments: 4 were expected but {len(sys.argv)} were given.")
        sys.exit(2)
    SAM_FILENAME = sys.argv[1]
    ALN_EXPL_OUTPUT_FILENAME = sys.argv[2]
    ERRORS_STATS_OUTPUT_FILENAME = sys.argv[3]
    NB_ALN_TO_COMPUTE = int(sys.argv[4])
    SPECIES_NAME = SAM_FILENAME.split("/")[-1].split(".sam")[0]

    # Used to further display progressing of the script
    STARTING_TIME = time.time()
    nb_aln_computed = 0

    ## Defines variables that will count sequencing errors
    # counters for number of match, mismatch and indel
    count_match, count_mismatch, count_insertion, count_deletion = 0, 0, 0, 0
    # dictionary that will store number of occurrences of each substitution
    substitution_occurrences_dict = initialize_substitution_dictionary()
    # dictionary that will store number of occurrences of insertions depending on their length
    insertion_length_occurrences_dict = {}
    # dictionary that will store number of occurrences of deletions depending on their length
    deletion_length_occurrences_dict = {}


    ## Computes explicit alignment and some basic sequencing error statistics
    ALN_EXPL_FILE = open(ALN_EXPL_OUTPUT_FILENAME, "w")
    SAM_FILE = open(SAM_FILENAME, "r")


    id_read = "" # initialize read ID

    for line in SAM_FILE:
        # skip useless first lines
        if line[0] == "@" or line == "\n":
            NB_ALN_TO_COMPUTE -= 1
            continue

        line_list = line.replace("\n", "").split("\t")


        # Display progression of script to STDOUT
        nb_aln_computed += 1
        prct_progressing = int(nb_aln_computed / NB_ALN_TO_COMPUTE * 100)
        if prct_progressing % 1 == 0:
            display_progressing_bar(prct_progressing, time.time())

        # check if the read is unmapped
        flag = line_list[1]
        if flag in ["4", "256", "2048"]:
            continue

        # check if this is not a secondary alignment
        read = line_list[9]

        if read == "0" or read == "*" or id_read == line_list[0]:
            continue

        # retrieve information about mapping
        id_read = line_list[0]
        if "reverse" in SAM_FILENAME:
            strand = " - (reverse)"
        else:
            strand = " + (forward)"
        position_in_genome = int(line_list[3])
        cigar_str = line_list[5]
        md_full_str = [elt for elt in line_list[-5:] if "MD:Z:" in elt][0]
        md_str = md_full_str.split("MD:Z:")[-1]



        # parse CIGAR, and update read and mapping position
        cigar, read, position_in_genome, soft_clips = parse_cigar(cigar_str, read, position_in_genome)



        # computes explicit alignment from read, CIGAR and MD
        alignment = Alignment(read, cigar, md_str)
        try:
            alignment.compute_explicit_alignment()
            alignment.checks_consistency()
        except Exception as exception:
            sys.stderr.write(f"The alignment for read {id_read} has been skipped because of\
                  following error: {exception}")
            continue


        # Update global dictionaries for insertions, deletions and substitutions
        substitution_occurrences_dict = update_dictionary(alignment.substitution_counts_dict,
                                                          substitution_occurrences_dict)
        insertion_length_occurrences_dict = update_dictionary(alignment.insertion_length_dict,
                                                              insertion_length_occurrences_dict)
        deletion_length_occurrences_dict = update_dictionary(alignment.deletion_length_dict,
                                                             deletion_length_occurrences_dict)

        # Update global counters of matches, mismatches, insertions and deletions
        count_match += alignment.match
        count_mismatch += alignment.mismatch
        count_insertion += alignment.insertion
        count_deletion += alignment.deletion

        # Write alignment in dedicated output file (header, genome and read)
        ALN_EXPL_FILE.write(" ; ".join([f"Read_ID={id_read}", f"strand={strand}",
                            f"position_in_genome={position_in_genome}",
                            f"bases_soft_clipped={soft_clips}"]) + "\n")
        ALN_EXPL_FILE.write(f"{alignment.genome_aln}\n{alignment.read_aln}\n")


    ALN_EXPL_FILE.close()
    SAM_FILE.close()

    # Write sequencing errors' basics statistics in dedicated output file
    ERRORS_STAT_FILE = open(ERRORS_STATS_OUTPUT_FILENAME, "w")

    sum_errors = sum([count_mismatch, count_insertion, count_deletion])
    total_alignment_length = sum([count_match, sum_errors])
    global_error_rate = round(100 - (count_match / total_alignment_length * 100), 4)
    mismatch_rate = round(count_mismatch / total_alignment_length * 100, 4)
    insertion_rate = round(count_insertion / total_alignment_length * 100, 4)
    deletion_rate = round(count_deletion / total_alignment_length * 100, 4)

    ERRORS_STAT_FILE.write(f"Total alignment length: {total_alignment_length}\n")
    ERRORS_STAT_FILE.write(f"Global error rate (%): {global_error_rate}\n")
    ERRORS_STAT_FILE.write(f"Insertion rate (%): {insertion_rate}\n")
    ERRORS_STAT_FILE.write(f"Deletion rate (%): {deletion_rate}\n")
    ERRORS_STAT_FILE.write(f"\nSubstitutions (genome to read):\nSubstitutions\tOccurrences\n")
    for substitution in substitution_occurrences_dict:
        ERRORS_STAT_FILE.write(f"{substitution}\t{substitution_occurrences_dict[substitution]}\n")
    ERRORS_STAT_FILE.write(f"\nInsertion lengths distributions:\n")
    ERRORS_STAT_FILE.write(f"Length\tOccurrences\n")
    for length in sorted(insertion_length_occurrences_dict):
        ERRORS_STAT_FILE.write(f"{length}\t{insertion_length_occurrences_dict[length]}\n")
    ERRORS_STAT_FILE.write(f"\nDeletion lengths distributions:\n")
    ERRORS_STAT_FILE.write(f"Length\tOccurrences\n")
    for length in sorted(deletion_length_occurrences_dict):
        ERRORS_STAT_FILE.write(f"{length}\t{deletion_length_occurrences_dict[length]}\n")

    ERRORS_STAT_FILE.close()

    total_elapsed_time = time.time() - STARTING_TIME
    total_elapsed_time = str(datetime.timedelta(seconds=round(total_elapsed_time)))
    sys.stdout.write(f"\nDone (elapsed time: {total_elapsed_time}).\n")
