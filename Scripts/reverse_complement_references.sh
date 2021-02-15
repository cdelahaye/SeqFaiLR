#!/bin/bash

echo ""
echo "====="
echo "Running $0"

# For each reference genome (in dedicated directory), computes the reverse complement, for further alignment
# Uses the fastx_toolkit

ref_dir="./Data/Reference_genomes" # input data = directory containing reference genomes

# Computes reverse complements
for ref_file in $ref_dir/*.fasta
do
  basename_ref_file=$(basename $ref_file)
  echo $ref_file

  # If it is an already reverse complemented file, skip it
  if [[ $basename_ref_file =~ "reverse" ]]; then
    continue
  fi

  # If a reverse complemented version of this file exists, skip it
  species_name=${basename_ref_file%.*}
  reverse_ref_file=${ref_dir}/${species_name}_reverse.fasta
  if [ -f $reverse_ref_file ]; then
    continue
  fi

  # Convert multiline fasta to one-line and convert to upper case
  temp_ref_file=${ref_dir}/${species_name}_temp.fasta
  python3 ./Scripts/Python/format_fasta.py $ref_file $temp_ref_file
  mv $temp_ref_file $ref_file

  # Computes the reverse complement and store it in new file
  paste -d "\n"  <(grep ">" $ref_file) <(grep -v ">" $ref_file | tr ATGCN TACGN | rev) > $reverse_ref_file

  echo "Done."
  echo ""

done

echo ""
echo "All reference genomes processed."
echo "---"
