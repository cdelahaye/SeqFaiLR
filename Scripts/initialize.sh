#!/bin/bash

#This script downloads needed extra algorithms
#  (minimap2, fastx_toolkit)
#  and run some prerequisites for further analyses

echo ""
echo "===== Running $0 ====="


ref_gen_dir="./Data/Reference_genomes" # input data: directory containing reference genomes
raw_read_dir="./Data/Raw_reads"

output_groups="./Data/species_groups.txt" # output file (groupID ; species) space separated
output_colors="./Data/species_color_GC.txt" # output file (species name ; color ; GC content of genome in %)


# --- --- ---
# Install minimap2 for alignments
echo "Check wether minimap2 is installed"
if ! command -v ./Scripts/minimap2/minimap2 &> /dev/null
then
  echo "  Installing minimap2 aligner:"
  git clone https://github.com/lh3/minimap2
  mv minimap2 ./Scripts/
  cd ./Scripts/minimap2 && make && cd ../..
  echo "    minimap2 installed."
else
  echo "  minimap2 already installed."
  echo "  version:"
  ./Scripts/minimap2/minimap2 --version
fi
if [ $? -eq 0 ]; then
  echo "Done."
else
  echo "Failed."
  exit 1
fi
echo ""


# --- --- ---
# For more readability, or to ease comparison, species data can be grouped, instead
#   of plotting all detailled results on graph.
# Put your raw read data in subfolder to do so
# If you do not want to group data, do not create subfolder and leave data as is in
#   ./Data/Raw_reads

# This script will create text file with information of groups, and will be used for
#  analyses

echo "Look at the organization of the species' folders"
if [ -f "$output_groups" ]; then
  echo "  File $output_groups already exist, do not compute it again (or remove it first)"
else
  python3 -u ./Scripts/Python/species_groups.py $raw_read_dir $output_groups
  if [ $? -eq 0 ]; then
    echo "Done."
  else
    echo "Failed."
    exit 1
  fi
fi
echo ""


# --- --- ---
# Compute GC content of reference genome for each species, and assign it a color
#  for future graphs

echo "Compute GC content of reference genome for each species, and assign it a color for future graphs"

if [ -f "$output_colors" ]; then
  echo "  File $output_colors already exists, do not compute it again (or remove it first)."
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
python3 -u ./Scripts/Python/species_color_GC.py $ref_gen_dir $output_colors

if [ $? -eq 0 ]; then
  echo "Done."
else
  echo "Failed."
  exit 1
fi
echo ""


# ----------
# Checks that every species data is associated to a reference genome
#  with identical naming
echo "Checks association species-genome and that files are correctly named"
python3 -u ./Scripts/Python/check_species_genome.py $ref_gen_dir $raw_read_dir
if [ $? -eq 0 ]; then
  echo "Done."
else
  echo "Failed."
  exit 1
fi
echo ""


# ----------
# For each reference genome, computes the reverse complement, for further alignment
# Uses the fastx_toolkit

echo "Computes reverse complements of each reference genome file"
# Computes reverse complements
for ref_file in $ref_gen_dir/*.fasta
do
  basename_ref_file=$(basename $ref_file)

  # If it's an already reverse complemented file, skip it
  if [[ $basename_ref_file =~ "reverse" ]]; then
    continue
  fi

  echo "  $basename_ref_file"

  # If a reverse complemented version of this file exists, skip it
  species_name=${basename_ref_file%.*}
  reverse_ref_file=${ref_gen_dir}/${species_name}_reverse.fasta
  if [ -f $reverse_ref_file ]; then
    continue
  fi

  # Convert multiline fasta to one-line and convert to upper case
  temp_ref_file=${ref_gen_dir}/${species_name}_temp.fasta
  python3 ./Scripts/Python/format_fasta.py $ref_file $temp_ref_file
  mv $temp_ref_file $ref_file

  # Computes the reverse complement and store it in new file
  paste -d "\n"  <(grep ">" $ref_file) <(grep -v ">" $ref_file | tr ATGCN TACGN | rev) > $reverse_ref_file

done

if [ $? -eq 0 ]; then
  echo "Done."
else
  echo "Failed."
  exit 1
fi

