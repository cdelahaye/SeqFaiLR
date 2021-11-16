#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Convert multi-line fasta to one line.
And also convert all characters to upper case.
"""

# --------------------------------------------------------------------------------------------------
# Packages

import sys

# --------------------------------------------------------------------------------------------------
# Main

if __name__ == "__main__":

    ## Parses arguments
    if len(sys.argv) != 3:
        print(f"ERROR: Wrong number of arguments: 2 expected but {len(sys.argv)-1} given.")
        sys.exit(2)
    REF_GEN_FILENAME = sys.argv[1]
    OUTPUT_FILENAME = sys.argv[2]

    file = open(REF_GEN_FILENAME, "r")
    new_file = open(OUTPUT_FILENAME, "w")

    first_line = file.readline()
    new_file.write(first_line)
    for line in file:
        if line[0] == ">":
            new_file.write("\n" + line)
            continue
        new_file.write(line.rstrip().upper())

    file.close()
    new_file.close()
