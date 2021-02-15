#!/bin/bash

echo ""
echo "====="
echo "Running $0"

# Compute GC content of reference genome for each species, and assign it a color for future graphs

ref_gen_dir="./Data/Reference_genomes" # input data: directory containing reference genomes
output="./Data/species_gc_color.txt" # output file (species name ; color ; gc content of genome)

if [ -f "$output" ]; then
  echo "$output already exists, skip this step."
  exit 0
fi

# Check that reference genome are given
if [ ! -d $ref_gen_dir ]; then
  echo "ERROR: Directory $ref_gen_dir cannot be found. Please provide reference genomes there." >&2
  exit 1
fi
if ! ls ${ref_gen_dir}/*.fasta &> /dev/null; then
  echo "ERROR: No fasta file of reference genome could be found in $ref_gen_dir. Please provide reference genomes in fasta format." >&2
  exit 1
fi

# Runs the script
python3 -u ./Scripts/Python/species_color_GC.py $ref_gen_dir $output
