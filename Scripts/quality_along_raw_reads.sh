#!/bin/bash

echo ""
echo "====="
echo "Running $0"

# From raw fastq files, plots quality along read (relative position).
# Also uses information from alignment to get mean number of bases clipped at read start and end

raw_reads_dir="./Data/Raw_reads/" # input data: directory containing fastq files of raw reads
aln_expl_dir="./Data/Alignment/" # input data: directory containing explicit alignment
output_graph_dir="./Output/Quality_along_reads/" # output directory: graphics
output_raw_dir="./Analysis/Quality_along_reads/" # output directory: raw results
file_color_gc_species="./Data/species_color_GC.txt" # file containing GC content of each species, and a dedicated color for plots
read_end_length=$1 # size of read ends (for plot)

mkdir -p $output_raw_dir
mkdir -p $output_graph_dir


# Check that fastq files and explicit alignment exists
if [ ! -d $raw_reads_dir ]; then
  echo "ERROR: Directory $raw_reads_dir cannot be found. Please provide fastq files in this directory." >&2
  exit 1
fi
if ! ls ${raw_reads_dir}/*.fastq &> /dev/null; then
  echo "ERROR: No fastq file could be found in $raw_reads_dir. Please provide fastq files in this directory." >&2
  exit 1
fi
if [ ! -d $aln_expl_dir ]; then
  echo "ERROR: Directory $aln_expl_dir cannot be found. Please run alignment_explicit.sh script prior to this one." >&2
  exit 1
fi
if ! ls ${aln_expl_dir}/*.txt &> /dev/null; then
  echo "ERROR: No explicit alignment could be found in $aln_expl_dir. Please run alignment_explicit.sh script prior to this one." >&2
  exit 1
fi


# Runs the script
python3 -u ./Scripts/Python/quality_along_raw_reads.py $raw_reads_dir $aln_expl_dir $output_raw_dir $output_graph_dir $file_color_gc_species $read_end_length

echo ""
echo "Done."
echo "---"
