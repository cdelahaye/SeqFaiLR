#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
For more readability, or to ease comparison, species data can be grouped, instead
  of plotting all detailled results on graph.
Put your raw read data in subfolder to do so
If you do not want to group data, do not create subfolder and leave data as is in ./Data/Raw_reads
If you want all your data to be merged on plots, create one subfolder and put all data into it.

This script will create text file with information of groups, and will be used for analyses
"""

# --------------------------------------------------------------------------------------------------
# Packages

import sys
import os

try:
    import seaborn as sns; sns.set_theme()
except ImportError:
    print("Package seaborn is required, please install it")
    sys.exit(1)

# --------------------------------------------------------------------------------------------------
# Main

if __name__ == "__main__":

    # Parse arguments
    NUMBER_EXPECTED_ARGUMENTS = 2
    if len(sys.argv) != NUMBER_EXPECTED_ARGUMENTS + 1:
        print(f"ERROR: Wrong number of arguments: {NUMBER_EXPECTED_ARGUMENTS} expected but "
              + f"{len(sys.argv)-1} given.")
        sys.exit(2)
    RAW_READ_DIR = sys.argv[1]
    OUTPUT_FILENAME = sys.argv[2]


    # Are path in RAW_READ_DIR files or subdirectories ?
    test_dir_or_file = ["dir" if os.path.isdir(RAW_READ_DIR + "/" + path)
                        else "file" if os.path.isfile(RAW_READ_DIR + "/" + path)
                        else "unknown"
                        for path in os.listdir(RAW_READ_DIR)]

    # Verify that there are only files or only subdirectories but not both
    if test_dir_or_file.count(test_dir_or_file[0]) != len(test_dir_or_file):
        print("  Found directories AND files. Please either put all your data into subdirectories, " +
              f"OR leave all your data in {RAW_READ_DIR} without making any subdirectory")
        sys.exit(1)


    # if not grouping of data -> leave it as is
    if test_dir_or_file[0] == "file":
        sys.exit(0)

    # if grouping data:
    elif test_dir_or_file[0] == "dir":
        nb_subdir = len(test_dir_or_file)
        list_colors = sns.diverging_palette(240, 10, n=nb_subdir, as_cmap=False)
        list_colors_hex = list_colors.as_hex() # get hexadecimal code for each color
        print(f"  Storing information in {OUTPUT_FILENAME}")
        output = open(OUTPUT_FILENAME, "w")
        for i, subdir in enumerate(os.listdir(RAW_READ_DIR)):
            for species in os.listdir(RAW_READ_DIR + "/" + subdir):
                output.write(f"{subdir}\t{species}\t{list_colors_hex[i]}\n")
        output.close()
