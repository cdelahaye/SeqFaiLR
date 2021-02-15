#!/bin/bash

#This script initiate structure for analysis
# and downloaded needed algorithms (for the moment
# minimap2 only is concerned).

echo ""
echo "====="
echo "Running $0"

# Creates folder to organise analysis
mkdir -p ./Data/Raw_reads
mkdir -p ./Data/Reference_genomes

# Install minimap2 for alignments
if [ ! -d ./Scripts/minimap2 ]; then
  echo "Installing minimap2 aligner:"
  git clone https://github.com/lh3/minimap2
  mv minimap2 ./Scripts/
  cd ./Scripts/minimap2 && make && cd ../
  echo "minimap2 installed."
fi

echo ""
echo "Initializing done. Please:"
echo "  - put raw reads data (in .fastq format*) in ./Data/Raw_reads"
echo "    (*: some analyse concern quality scores, thus can only be performed with .fastq files)"
echo "  - put reference genomes in ./Data/Reference_genomes"
echo "  WARNING: raw reads and reference genomes are expected to be named similarly, for example escherichia_coli.fastq and escherichia_coli.fasta respectively"
echo ""
