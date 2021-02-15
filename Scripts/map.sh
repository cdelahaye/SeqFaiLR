#!/bin/bash

echo ""
echo "====="
echo "Running $0"

# This script uses minimap2 aligner to map raw reads against their associated reference genome
# The following line ensures minimap2 is present in the directory
if [ ! -d ./Scripts/minimap2 ]; then
  echo "ERROR: minimap2 cannot be found. Please run ./initialize.sh script." >&2
  exit 1
fi


# 1/ Mapping reads
echo "Step 1: Mapping reads"
output_path=./Data/Mapping

for path_raw_read_file in ./Data/Raw_reads/*.fastq
do
  raw_read_file=$(basename $path_raw_read_file)
  species_name=${raw_read_file%.*}

  # Try to map on forward strand of reference genome
  reference_genome_file=./Data/Reference_genomes/${species_name}.fasta
  if [ ! -f $reference_genome_file ]; then
    echo "  ERROR (file not found). Expected to find $reference_genome_file for $species_name raw read but such file does not exist. Please check spelling of this file." >&2
    exit 1
  fi
  output_file=${output_path}/${species_name}.sam
  mkdir -p $output_path
  if [ -f $output_file ]; then
    echo "SKIPPED ALIGNMENT"
    echo "Alignement file already exists for species ${species_name}, so alignment will not be performed again. If you want the alignment to be performed, please remove $output_file first, then re-run this script."
    echo ""
    continue
  fi
  echo "Mapping reads for $species_name ($path_raw_read_file vs $reference_genome_file), on forward strand"
  echo "Output .sam file will be stored in $output_file"
  ./Scripts/minimap2/minimap2 --MD --eqx -ax map-ont --for-only --secondary=no --sam-hit-only $reference_genome_file $path_raw_read_file > $output_file
  echo "Done."
  echo ""

  # Try to map on reverse strand of reference genome
  reference_genome_file=./Data/Reference_genomes/${species_name}_reverse.fasta
  if [ ! -f $reference_genome_file ]; then
    echo "  ERROR (file not found). Expected to find $reference_genome_file for $species_name raw read but such file does not exist. Please check spelling of this file, and check that the reference genome's reverse complement has been computed." >&2
    exit 1
  fi
  output_file=${output_path}/${species_name}_reverse.sam
  mkdir -p $output_path
  if [ -f $output_file ]; then
    echo "SKIPPED ALIGNMENT"
    echo "Alignement file already exists for species ${species_name}, so alignment will not be performed again. If you want the alignment to be performed, please remove $output_file first, then re-run this script."
    echo ""
    continue
  fi
  echo "Mapping reads for $species_name ($path_raw_read_file vs $reference_genome_file)"
  echo "Output .sam file will be stored in $output_file"
  ./Scripts/minimap2/minimap2 --MD --eqx -ax map-ont --for-only --secondary=no --sam-hit-only $reference_genome_file $path_raw_read_file > $output_file
  echo "Done."
  echo ""

done


# 2/ Cleaning alignment sam files
echo ""
echo "Step 2: Cleaning sam file for wrong mappings"
echo ""

sam_dir=./Data/Mapping
for sam_file in $sam_dir/*_reverse.sam
do
  sam_file=${sam_file/_reverse/}
  basename_sam_file=$(basename $sam_file)
  species_name=${basename_sam_file%.*}
  echo $species_name
  python3 -u ./Scripts/Python/clean_sam.py $species_name $sam_dir
done


# Replace old alignment files by their cleaned versions
sam_dir=./Data/Mapping
for sam_file_reverse_clean in $sam_dir/*_reverse_clean.sam
do
  echo $sam_file_reverse_clean
  sam_file_forward_clean=${sam_file_reverse_clean/_reverse/}
  sam_file=${sam_file_forward_clean/_clean/}
  basename_sam_file=$(basename $sam_file)
  species_name=${basename_sam_file%.*}
  echo $species_name
  sam_file_forward=${sam_file_forward_clean/_clean/}
  mv $sam_file_forward_clean $sam_file_forward
  sam_file_reverse=${sam_file_reverse_clean/_clean/}
  mv $sam_file_reverse_clean $sam_file_reverse
done
