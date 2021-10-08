#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
For more readability, or to ease comparison, species data can be grouped, instead
  of plotting all detailled results on graph.
It will ask you a name (character or numbers, please avoid special character and space)
If you do not want to group data, answer "NO"
If you want to group all data together, answer "ALL"
You can still edit the file afterwards.

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

    dict_group_species = {}

    # Ask if want to group all or group nothing
    print("Do you want to group species on plots?")
    print("  - If you want to group all species together, type ALL")
    print("  - If you do not want to group species, type NO")
    print("  - If you want to create multiple groups, type anything else!")
    do_group_all_none = input("").upper()

    if do_group_all_none == "ALL":
        print(" Okay, I'll merge everything!")
        group_name = "All data"
        dict_group_species[group_name] = []
    elif do_group_all_none == "NO":
        print(" Okay, I wont merge anything!")
        nb_color = len(os.listdir(RAW_READ_DIR))
    else:
        print(" Okay, now tell me more about it!")
    print("")


    for filename in os.listdir(RAW_READ_DIR):
        if do_group_all_none == "ALL":
            dict_group_species[group_name] += [filename]
        elif do_group_all_none == "NO":
            dict_group_species[filename] = [filename]
        else:
            group_name = input(f"  Which group for {filename}?\t")
            if group_name not in dict_group_species:
                dict_group_species[group_name] = []
            dict_group_species[group_name] += [filename]

    nb_color = len(dict_group_species)
    list_colors = sns.diverging_palette(240, 10, n=nb_color, as_cmap=False)
    list_colors_hex = list_colors.as_hex() # get hexadecimal code for each color

    output = open(OUTPUT_FILENAME, "w")

    for i,group in enumerate(dict_group_species):
        color = list_colors_hex[i]
        for species in dict_group_species[group]:
            output.write(f"{group}\t{species}\t{color}\n")
    output.close()
