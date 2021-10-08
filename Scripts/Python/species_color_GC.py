#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Compute GC content of reference genome for each species.
Then associate each species with a color (for graphs) and gc content.
Output all of this in tab separated file
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
        print(f"ERROR: Wrong number of arguments: {NUMBER_EXPECTED_ARGUMENTS} expected but {len(sys.argv)-1} given.")
        sys.exit(2)
    REF_DIRNAME = sys.argv[1]
    OUTPUT = sys.argv[2]

    if REF_DIRNAME[-1] != "/":
        REF_DIRNAME += "/"

    dict_gc_species = {}

    for ref_filename in os.listdir(REF_DIRNAME):

        if "_reverse" in ref_filename:
            continue # For other scripts' needs, reverse complement of reference genome may have already been computed ; ignore them

        species_name = ref_filename.split(".")[0]
        gc_count, total_bases_count = 0, 0
        file = open(REF_DIRNAME + ref_filename, "r")
        for line in file:
            if line[0] == ">":
                continue
            line = line.rstrip().upper().replace("N", "") # Ignore N bases
            gc_count += line.count("G") + line.count("C")
            total_bases_count += len(line)
        file.close()

        gc_content = gc_count / total_bases_count * 100
        if gc_content not in dict_gc_species:
            dict_gc_species[gc_content] = []
        dict_gc_species[gc_content] += [species_name]


    nb_species = len(dict_gc_species.values())

    list_colors = sns.diverging_palette(240, 10, n=nb_species, as_cmap=False)
    list_colors_hex = list_colors.as_hex() # get hexadecimal code for each color
    #sns.palplot(list_colors) # to visualise

    file = open(OUTPUT, "w")
    for i, gc_content in enumerate(sorted(dict_gc_species.keys())):
        color = list_colors_hex[i]
        for species_name in dict_gc_species[gc_content]:
            file.write(f"{species_name}\t{color}\t{round(gc_content, 2)}\n")

    file.close()
