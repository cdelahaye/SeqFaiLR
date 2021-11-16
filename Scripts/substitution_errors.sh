#!/bin/bash

echo ""
echo "====="
echo "Running $0"

# From outputs of errors, plots relative abundance of each type of substitution

errors_dir="./Analysis/Error_rates" # input data: directory containing substitution errors counts
output_dir="./Output/"
output_graph="./Output/substitution_errors.png" # output filename
path_file_color_gc_species="./Data/species_color_GC.txt" # file containing GC content of each species
path_file_species_groups="./Data/species_groups.txt" # file containing species groups

mkdir -p $output_dir

# Check that input files exist
if [ ! -d $errors_dir ]; then
  echo "ERROR: Directory $error_dir cannot be found. Please run alignment_explicit.sh script prior to this one." >&2
  exit 1
fi
if ! ls ${errors_dir}/*.txt &> /dev/null; then
  echo "ERROR: No file could be found in $error_dir. Please run alignment_explicit.sh script prior to this one." >&2
  exit 1
fi


# Runs the script
python3 -u ./Scripts/Python/substitution_errors.py $errors_dir $output_graph $path_file_color_gc_species $path_file_species_groups

echo "Done."
echo "---"
