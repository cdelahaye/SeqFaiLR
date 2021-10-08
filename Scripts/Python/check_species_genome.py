#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Verifies that each species has an associated reference genome, and vice versa.
"""

# --------------------------------------------------------------------------------------------------
# Packages

import sys
import os


# --------------------------------------------------------------------------------------------------
# Main

if __name__ == "__main__":

    # Parse arguments
    NUMBER_EXPECTED_ARGUMENTS = 2
    if len(sys.argv) != NUMBER_EXPECTED_ARGUMENTS + 1:
        print(f"ERROR: Wrong number of arguments: {NUMBER_EXPECTED_ARGUMENTS} expected but "
              + f"{len(sys.argv)-1} given.")
        sys.exit(2)
    REF_GENOME_DIR = sys.argv[1]
    RAW_READS_DIR = sys.argv[2]


    dict_genome_species = {}

    for ref_genome in os.listdir(REF_GENOME_DIR):
        ref_genome = ref_genome.replace(".fasta", "").replace("_reverse", "")
        if ref_genome not in dict_genome_species:
            dict_genome_species[ref_genome] = []

    for path in os.listdir(RAW_READS_DIR):
        if os.path.isdir(RAW_READS_DIR + "/" + path):
            for filename in os.listdir(RAW_READS_DIR + "/" + path):
                species_name = filename.split(".f")[0] # work for .fasta .fastq .fa .fq ...
                if species_name in dict_genome_species:
                    dict_genome_species[species_name] += [filename]
                else:
                    full_path = RAW_READS_DIR + "/" + path + filename
                    print(f"ERROR: no reference genome found for species {full_path}",
                          "Please check spelling of files so that XXX.fasta reference genome exactly",
                          "matches with XXX.fastq species name")
                    sys.exit(1)
        else:
            species_name = path.split(".f")[0]
            if species_name in dict_genome_species:
                dict_genome_species[species_name] += [filename]
            else:
                full_path = RAW_READS_DIR + "/" + path
                print(f"ERROR: no reference genome found for species {full_path}",
                      "Please check spelling of files so that XXX.fasta reference genome exactly",
                      "matches with XXX.fastq species name")
                sys.exit(1)

    for ref_genome in dict_genome_species:
        if dict_genome_species[ref_genome] == []:
            print(f"ERROR: no species associated to reference genome {ref_genome}.",
                      "Please check spelling of files so that XXX.fasta reference genome exactly",
                      "matches with XXX.fastq species name, or remove this reference genome if unused.")
            sys.exit(1)
