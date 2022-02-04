#!/bin/bash

echo ""
echo "====="
echo "Running $0"

# From SAM files, computes an explicit alignment, i.e. provides the sequences of aligned genome part and read"

# This program retrieves alignment of reads against a genome, by using MD field and CIGAR string provided in the .sam file.
# It outputs the final alignments with, for each read, its alignment against the genome (the genome is first printed, then the read)"

# Output format :
#  - (1) Header that contains information about the read aligned
#  - (2) the alignment of the genome part
#  - (3) the alignment of the read"

# It also computes substitution, insertion and deletion errors ; which are stored in a file called *_err_substi_indel.txt."

map_dir="./Data/Mapping/" # input data
aln_expl_dir="./Data/Alignment/" # output data containing explicit alignment
seq_err_stat_dir="./Analysis/Error_rates/" # output data containing basic information about sequencing error rates

# Check that mapping data (.sam files) exists
if [ ! -d $map_dir ]; then
  echo "ERROR: Directory $map_dir cannot be found. Please run map.sh before this script." >&2
  exit 1
fi
if ! ls ${map_dir}/*.sam &> /dev/null; then
  echo "ERROR: No SAM file could be found in $map_dir. Please run map.sh before this script." >&2
  exit 1
fi

# Creates output directory that will contain explicit alignments
if [ ! -d $aln_expl_dir ]; then
  mkdir -p $aln_expl_dir
fi

# Creates output directory that will contain sequencing error basic information
if [ ! -d $seq_err_stat_dir ]; then
  mkdir -p $seq_err_stat_dir
fi




# Computes alignment
echo "Compute alignments"
for sam_file in $map_dir/*.sam
do
  basename_sam_file=$(basename $sam_file)
  species_name=${basename_sam_file%.*}
  aln_expl_file=${aln_expl_dir}/${species_name}.txt
  seq_err_stat_file=${seq_err_stat_dir}/${species_name}.txt
  if [ -f ${aln_expl_file/_reverse/} ]; then
    echo "Explicit alignement file already exists for species ${species_name}, so it will not be computed again. If you want the alignment to be performed, please remove $aln_expl_file first, then re-run this script"
    echo ""
    continue
  fi
  echo -n "  - $species_name"
  nb_aln_to_compute=$(wc -l < $sam_file)
  python3 ./Scripts/Python/alignment_explicit.py $sam_file $aln_expl_file $seq_err_stat_file $nb_aln_to_compute
done

if [ $? -eq 0 ]; then
  echo "Done."
else
  echo "Failed."
  exit 1
fi

echo "---"


echo ""
echo "Merging of forward and reverse alignment file for each species"

# checks if some *reverse files exist: if true, then they have to be merged
if ls $aln_expl_dir/*reverse.txt > /dev/null 2>&1; then
  for aln_file_reverse in ${aln_expl_dir}*reverse.txt
  do
    echo "  - $aln_file_reverse"
    basename_aln_file=$(basename $aln_file_reverse)
    species_name=${basename_aln_file%%"_reverse"*}
    aln_file_forward=${aln_expl_dir}/${species_name}.txt
    tmp_file=${aln_expl_dir}/${species_name}
    cat $aln_file_reverse $aln_file_forward > $tmp_file
    mv $tmp_file $aln_file_forward
    rm $aln_file_reverse
  done
else
  echo "  Already done, nothing to merge."
fi


echo ""
echo "---"

echo "Merging error stat files (merging of forward and reverse results)"

if ls $seq_err_stat_dir/*reverse* > /dev/null 2>&1; then
  for stat_err_file_reverse in ${seq_err_stat_dir}*reverse*
  do
    echo $stat_err_file_reverse
    stat_err_file_forward=${stat_err_file_reverse/_reverse/}
    temp_file=${stat_err_file_forward/.txt/_temp.txt}
    python3 -u ./Scripts/Python/merge_error_stat.py $stat_err_file_forward $stat_err_file_reverse $temp_file
    rm $stat_err_file_reverse
    mv $temp_file $stat_err_file_forward
  done
else
  echo "  Already done, nothing to merge."
fi


