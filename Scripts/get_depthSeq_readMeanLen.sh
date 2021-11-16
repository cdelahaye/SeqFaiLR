#!/bin/bash

echo ""
echo "====="
echo "Running $0"

# From alignment files, computes: reference genome length, sequencing depth, GC content in reference genome, mean length of sequenced reads
# Output result in raw file, both in tabulated and LaTeX format

aln_expl_dir="./Data/Alignment/" # input data: directory containing explicit alignment
ref_dir="./Data/Reference_genomes/" # input data: directory of reference genomes
output_raw_dir="./Analysis/Sequencing_depth_read_length/" # output directory: raw results

mkdir -p $output_raw_dir

# Check that fastq files and explicit alignment exists
if [ ! -d $aln_expl_dir ]; then
  echo "ERROR: Directory $aln_expl_dir cannot be found. Please run alignment_explicit.sh script prior to this one." >&2
  exit 1
fi
if ! ls ${aln_expl_dir}/*.txt &> /dev/null; then
  echo "ERROR: No explicit alignment could be found in $aln_expl_dir. Please run alignment_explicit.sh script prior to this one." >&2
  exit 1
fi


# Runs the script
python3 -u ./Scripts/Python/get_depthSeq_readMeanLen.py $aln_expl_dir $ref_dir $output_raw_dir

echo ""
echo "Done."
echo "---"
