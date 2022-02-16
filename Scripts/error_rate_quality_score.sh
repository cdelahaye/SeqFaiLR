#!/bin/bash

echo ""
echo "====="
echo "Running $0"

# From alignment files, plots the error rate of reads depending on their quality scores

aln_expl_dir="./Data/Alignment/" # input data: directory containing explicit alignment
raw_read_dir="./Data/Raw_reads/" # input data: directory containing fastq files
output_graph_dir="./Output/Error_rate_quality_score/" # output directory: graphics
output_raw_dir="./Analysis/Error_rate_quality_score/" # output directory where raw results are stored
file_color_gc_species="./Data/species_color_GC.txt" # file containing GC content of each species, and a dedicated color for plots
file_species_groups="./Data/species_groups.txt" # file containing species groups

window_length=$1
err_qual_max_aln_nb=$2
err_qual_min_occ_read=$3
err_qual_min_occ_window=$4


mkdir -p $output_graph_dir
mkdir -p $output_raw_dir

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
python3 -u ./Scripts/Python/error_rate_quality_score.py $aln_expl_dir $raw_read_dir $output_raw_dir $output_graph_dir $file_color_gc_species $file_species_groups $window_length $err_qual_max_aln_nb $err_qual_min_occ_read $err_qual_min_occ_window

echo "Done."
echo "---"
