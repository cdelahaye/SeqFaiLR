#!/bin/bash

echo ""
echo "====="
echo "Running $0"

# Analyse explicit alignment files to compute error rates depending on GC content of reads.

aln_expl_dir="./Data/Alignment/" # input data: directory containing explicit alignment
output_graph_dir="./Output/Error_rate_gc_read/" # output directory: graphics
output_raw_dir="./Analysis/Error_rate_gc_read/" # output directory: raw results
file_color_gc_species="./Data/species_color_GC.txt" # file containing GC content of each species, and a dedicated color for plots
file_species_groups="./Data/species_groups.txt" # file containing species groups
nb_max_aln=$1
min_occ=$2

mkdir -p $output_raw_dir
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
python3 -u ./Scripts/Python/error_rate_reads_gc_content.py $aln_expl_dir $output_raw_dir $output_graph_dir $file_color_gc_species $file_species_groups $nb_max_aln $min_occ

echo "Done."
echo "---"
