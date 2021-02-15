#!/bin/bash

echo ""
echo "====="
echo "Running $0"

# From alignment files, plots the relative coverage depending on GC content

aln_expl_dir="./Data/Alignment" # input data: directory containing explicit alignment
ref_dir="./Data/Reference_genomes" # input data: directory containing reference genomes
output_graph_dir="./Output/Coverage_gc_bias" # output directory: graphics
path_file_color_gc_species="./Data/species_gc_color.txt" # file containing GC content of each species, and a dedicated color for plots

mkdir -p $output_graph_dir


# Check that explicit alignment exists
if [ ! -d $aln_expl_dir ]; then
  echo "ERROR: Directory $aln_expl_dir cannot be found. Please run alignment_explicit.sh script prior to this one." >&2
  exit 1
fi
if ! ls ${aln_expl_dir}/*.txt &> /dev/null; then
  echo "ERROR: No explicit alignment could be found in $aln_expl_dir. Please run alignment_explicit.sh script prior to this one." >&2
  exit 1
fi


# Runs the script
python3 -u ./Scripts/Python/coverage_GC_content.py $aln_expl_dir $ref_dir $path_file_color_gc_species $output_graph_dir

echo ""
echo "Done."
echo "---"
