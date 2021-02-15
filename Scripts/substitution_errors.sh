#!/bin/bash

echo ""
echo "====="
echo "Running $0"

# From outputs of errors, plot relative abundance of each type of substitution

errors_dir="./Analysis/Error_rates" # input data: directory containing substitution errors counts
output_graph_dir="./Output/Substitution_errors" # output directory: graphics
path_file_color_gc_species="./Data/species_gc_color.txt" # file containing GC content of each species, and a dedicated color for plots

mkdir -p $output_graph_dir

# Check that input files
if [ ! -d $errors_dir ]; then
  echo "ERROR: Directory $error_dir cannot be found. Please run alignment_explicit.sh script prior to this one." >&2
  exit 1
fi
if ! ls ${errors_dir}/*.txt &> /dev/null; then
  echo "ERROR: No file could be found in $error_dir. Please run alignment_explicit.sh script prior to this one." >&2
  exit 1
fi


# Runs the script
python3 -u ./Scripts/Python/substitution_errors.py $errors_dir $output_graph_dir $path_file_color_gc_species

echo ""
echo "Done."
echo "---"
