#!/bin/bash

echo ""
echo "====="
echo "Running $0"

# From explicit alignment files and reference genome files, computes a graph of error rates (detailed, and global)
#   depending on relative position in the reference genome.

sam_dir="./Data/Mapping/" # input data: directory containing explicit alignment
ref_gen_dir="./Data/Reference_genomes/" # input data: directory containing reference genomes
output_graph_dir="./Output/Error_rates_along_genome/" # output directory: graphics
file_color_gc_species="./Data/species_color_GC.txt" # file containing GC content of each species

mkdir -p $output_graph_dir


# Check that reference genomes and alignments exist
if [ ! -d $ref_gen_dir ]; then
  echo "ERROR: Directory $ref_gen_dir cannot be found. Please provide reference genomes there." >&2
  exit 1
fi
if ! ls ${ref_gen_dir}/*.fasta &> /dev/null; then
  echo "ERROR: No fasta file of reference genome could be found in $ref_gen_dir. Please provide reference genomes in fasta format." >&2
  exit 1
fi
if [ ! -d $sam_dir ]; then
  echo "ERROR: Directory $sam_dir cannot be found. Please run map.sh script prior to this one." >&2
  exit 1
fi
if ! ls ${sam_dir}/*.sam &> /dev/null; then
  echo "ERROR: No sam file could be found in $saml_dir. Please run map.sh script prior to this one." >&2
  exit 1
fi


# Runs the script
python3 -u ./Scripts/Python/error_rates_along_genome.py $sam_dir $ref_gen_dir $output_graph_dir $file_color_gc_species

echo ""
echo "Done."
echo "---"
