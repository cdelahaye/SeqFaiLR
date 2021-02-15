#!/bin/bash

# Parse command-line options

function show_usage (){
    echo "Usage $0 [options]"
    echo ""
    echo "Options:"
    echo "-h|--help, display this help"
    echo "-m|--map, align reads against their reference genome"
    echo "-l|--low_complexity, analyses errors on low complextiy regions: hoompolymers, heteropolymers and trinucleotides"
    echo "   --homopolymers, analyses errors on homopolymeric regions"
    echo "   --heteropolymers, analyses errors for heteropolymer patterns"
    echo "   --trinucleotides, analyses errors for trinucleotides patterns"
    echo "-q|--quality, computes quality scores along reads and errors rates depending on quality score"
    echo "-g|--gc, computes coverage and error rates of reads depending on their GC content"
    echo "-a|--all, run all these scripts"

return 0
}


if [[ $# -eq 0 ]];then
   show_usage
   exit 1
fi
if [[ $# -gt 1 ]];then
   echo "WARNING: only one argument expected, so only the first one will be processed"
   echo "If you want to run run all analysis, please use '$0 --all' or '$0 -a'"
   echo ""
fi


case "$1" in
     -h | --help )
       show_usage
       ;;
     -m | --map )
       mode=mapping
       ;;
     -l | --low_complexity )
       mode=low_complexity
       ;;
     --homopolymers )
       mode=homopolymers
       ;;
     --heteropolymers )
       mode=heteropolymers
       ;;
     --trinucleotides )
       mode=trinucleotides
       ;;
     -q | --quality )
       mode=quality
       ;;
     -g | --gc )
       mode=gc
       ;;
     -a | --all )
       mode=all
       ;;
     *)
      echo "Incorrect input provided"
      show_usage
      ;;
esac



# Some scripts are written in Python, the following lines ensures python is sourced, and of compatible version
python_version=$(python3 -V 2>&1 | grep -Po '(?<=Python )(.+)')
if [[ -z "$python_version" ]]
then
    echo "ERROR: No version of python could be found." >&2
    exit 1
fi
python_parsed_version=$(echo "${python_version//./}")
if ! [[ "$python_parsed_version" -gt "360" ]]
then
    echo "Error: Python version > 3.6 required"
    exit 1
fi






# Run appropriate scripts

case $mode in
   mapping )
     ./Scripts/initialize.sh
     ./Scripts/species_color_GC.sh
     ./Scripts/reverse_complement_references.sh
     ./Scripts/map.sh
     ./Scripts/alignment_explicit.sh
     ./Scripts/substitution_errors.sh
     ./Scripts/error_rates_along_genome.sh
     ./Scripts/get_depthSeq_readMeanLen.sh
     ;;
   quality )
     ./Scripts/error_rate_quality_score.sh
     ./Scripts/quality_along_raw_reads.sh
     ;;
   gc )
     ./Scripts/coverage_gc_content.sh
     ./Scripts/error_rates_gc_read.sh
     ;;
   low_complexity )
     ./Scripts/homopolymer_errors.sh
     ./Scripts/heteropolymer_errors.sh
     ./Scripts/trinucleotides.sh
     ;;
   homopolymers )
     ./Scripts/homopolymer_errors.sh
     ;;
   heteropolymers )
     ./Scripts/heteropolymer_errors.sh
     ;;
   trinucleotides )
     ./Scripts/trinucleotides.sh
     ;;
   all )
     ./Scripts/initialize.sh
     ./Scripts/species_color_GC.sh
     ./Scripts/reverse_complement_references.sh
     ./Scripts/map.sh
     ./Scripts/alignment_explicit.sh
     ./Scripts/substitution_errors.sh
     ./Scripts/error_rates_along_genome.sh
     ./Scripts/get_depthSeq_readMeanLen.sh
     ./Scripts/error_rate_quality_score.sh
     ./Scripts/quality_along_raw_reads.sh
     ./Scripts/coverage_gc_content.sh
     ./Scripts/error_rates_gc_read.sh
     ./Scripts/homopolymer_errors.sh
     ./Scripts/heteropolymer_errors.sh
     ./Scripts/trinucleotides.sh
     ;;
   * )
     echo "Not recognized"
     exit 1
     ;;
esac
