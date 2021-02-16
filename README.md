Algorithms used to analyse Nanopore long reads sequencing error profile.
The algorithms have been designed for Nanopore data, but can be applied for other long read data (and short read data too, however the analysis may not be relevant).

 
From raw reads and reference genomes, these scripts performs alignment and compute several analysis (low-complexity regions sequencing accuracy, GC bias, links between error rates and quality scores, and so one).

---

### Python requirements

Python > 3.6
numpy, matplotlib, os, sys, re, time, datetime, ast, seaborn


---

# How to run it

Here is describe how to run the pipeline. For details on each scripts, see below.

To display help: `seqfailr --help` or `seqfailr -h`.

### 1/ Aligning reads on their associated reference genome: `seqfailr --map` 

Short option: `seqfailr -m`

This command will run multiple scripts in order to:
  - download minimap2
  - assign a color to each species for future analysis and graphs
  - mapping reads against their reference genome
  - analyse substitution errors
  - compute error rate along reference genome (only available for species whose reference genome is one chromosome)
  - compute depth of coverage and mean read length

### 2/ Quality analysis: `seqfailr --quality`

Short option: `seqfailr -q`

This command will run scripts in order to compute:
  - link between error rate and quality scores (compared to expected Phred score)
  - quality scores along raw reads


### 3/ GC bias analysis: `seqfailr --gc`

Short option: `seqfailr -g`

This command will run scripts in order to compute:
  - coverage depending on GC content
  - error rates depending on GC content



### 4/ Analysis of low-complexity sequences `seqfailr --low_complexity`

Short option: `seqfailr -l`

Compute analysis for homopolymers ( (X)^n ), heteropolymers ( (XY)^n ) and trinucleotides ( (XYZ)^n ).
For individual run of these analysis, use:
  - `seqfailr --homopolymers`, results stored in Output/Homopolymer_errors and Analysis/Homopolymer_errors
  - `seqfailr --heteropolymers`, results stored in Output/Heteropolymer_errors and Analysis/Heteropolymer_errors
  - `seqfailr --trinucleotides`, results stored in Output/Trinucleotides and Analysis/Trinucleotides



User can use `seqfailr --all` or `seqfailr -a` to run all of this analysis at once.

---

# Toy example

Some (part of) fastq files, and their associate reference genome, are available in Data/Raw_reads and Data/Reference_genomes to try this pipeline.


---

# Details of scripts

The following describes main scripts.


### Basecalling

The basecaller used in the article was the ONT Guppy basecaller, version 4.2.2.
This basecaller is proprietary thus it cannot be provided here (a customer account in Nanopore is required).
Users have to run the basecalling step by themselves.



### (1) Mapping reads on genomes `map.sh`

For aligning Nanopore long reads we chose to use minimap2 aligner. 
A mode is dedicated to nanopore long read alignment (`-ax map-ont`). If you run our pipeline on other (long or short reads) data, please refer to the Minimap2 manpage (https://lh3.github.io/minimap2/minimap2.html) to choose appropriate settings.
In order to be able to reconstruct an explicit alignment (see next section) we also used parameters to retrieve MD tag (`--MD`) and CIGAR operators (`--eqx`) in SAM format (`-a`).

Minimap2 reference article: "Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018 Sep 15;34(18):3094-3100. doi: 10.1093/bioinformatics/bty191. PMID: 29750242; PMCID: PMC6137996."

#### Important note on dealing with reads mapping on reverse strand:
By default, if one read maps on the reverse strand of the reference genome (i.e. the reverse complement of the given genome in fasta file), the aligner maps the reverse complement of the read on the forward strand.
However, this raises two problems (we note R the sequenced read, R_rev its reverse complement, G the given reference genome and G_rev its reverse complement):
  - given that the aligner uses heuristics, computing alignment of R against G_rev, and computing alignment of R_rev against G, are not identical
  - even if such alignments were identical, it is important (for our study of sequencing errors) to keep exactly the information of the sequenced read: taking its reverse complement does not represent the reality of the sequencing
Thus, we try explicitely to align reads on the forward, and reverse complement of the reference genome, and we prevent the aligner to reverse complement the reads.
To do so, we compute the reverse complement of each reference genome, and store them in dedicated files (TODO: nom du script).


The command used to run minimap2 is the following:  
`
minimap2 --MD --eqx -ax map-ont --for-only --secondary=no --sam-hit-only $reference_genome_file $raw_reads_file
`
See the manpage of minimap2 for explanation of options: https://lh3.github.io/minimap2/minimap2.html

Mapping reads is handled by script `map.sh`. It will browse raw reads files (in Data/Raw_reads directory), will look for associated reference genomes (Data/Reference_genomes directory) and will output SAM alignment file (Data/Map directory).


### (1) Compute explicit alignment `alignment_explicit.sh`

From SAM files, computes the explicit alignment between reads and genome. 
/!\ MD and CIGAR tags in SAM file are mandatory for this script to run /!\  

The `alignment_explicit.sh` script browses each files in Data/Map and computes the alignment that is stored in Data/Alignment (skip if a corresponding file already exists in this directory). 
It also computes some basic statistics: total aligned bases for each file, error rates (global, and detailled for mismatches, insertions and deletions), occurrences of each substitution type, and distribution length of insertions and of deletions.


### (1) Compute ratio of substitution errors `substitution_errors.sh`

From results of explicit alignment, computes and plot ratio of all possible substitution errors.
Result graph is stored in /Output/Substitution_errors.


### (1) Compute error rates along genome `error_rates_along_genome.sh`

Only possible for species whose reference genome is in one part.
From alignment results, plot error rates depending on relative position in the reference genome.

Outputs are saved in /Output/Error_rates_along_genome (for graph) and /Analysis/Error_rates_along_genome (for raw results)

### (1) Compute depth of coverage and mean read length `get_depthSeq_readMeanLen.sh`

For each species, compute depth of coverage and its mean read length. Output is stored in Analysis/Sequencing_depth_read_length.


### (2) Error rate depending on quality scores `error_rate_quality_score.sh`

From raw reads extracts quality scores (fastq files are thus mandatory), and from explicit alignment extracts error rates.
Compute data to plot quality scores depending on error rate, both at read scale and at local scale (default, 100-bases windows). 

Plots are stored in Output/Error_rate_quality_score, and raw resuts (txt files) are stored in Analysis/Error_rate_quality_score.

### (2) Quality scores along raw reads `quality_along_raw_reads.sh`

From fastq files, compute and plot quality scores depending on the position considered in the raw read.
Plots are stored in Output/Quality_along_reads and raw results in Analysis/Quality_along_reads.


### (3) Depth of coverage depending on GC content `coverage_gc_content.sh`

From alignment files, plots the relative coverage depending on GC content.
Plots are stored in Output/Coverage_gc_bias.

### (3) Error rates depending on GC content `error_rates_gc_read.sh`

Analyse explicit alignment files to compute error rates depending on GC content of reads.
Plots are stored in Output/Error_rate_gc_read, and raw results in Analysis/Error_rate_gc_read.



### (4) Errors associated to homopolymer regions `homopolymer_errors.sh`

From reference genomes, computes distribution of homopolymers depending on the base of homopolymer (A, C, G or T).
Computes a plot (stored in Output/Homopolymer_errors/homopolymer_genomic_distribution.png) and raw results containing raw tab-separated table and a LaTeX template of these data (stored in Analysis/Homopolymer_errors/homopolymer_genomic_distribution.txt).

From explicit alignment files, computes:
 - global quantification of homopolymer errors (quantify_errors_homopolymer_vs_total.txt, no graph here), i.e. number of mismatches and indels in homopolymeric regions (in reference genome) and in all regions of the genome.
 - error rates (mismatches and indels) depending on homopolymer length (error_rates_depending_on_homopolymer_length*)
 - ratio of correctly sequenced homopolymers over total number of homopolymers, depending on their length (percentage_homopolymer_correctly_sequenced_by_length*). Results are presented with a boxplot where each box represent variations between species.
 - differences between expected homopolymer length (in reference genome) and sequenced one (difference_expected_sequenced_homopolymer_length*). Results are displayed as boxplots for which each box merges results for all species analysed.

### (4) Errors associated to heteropolymer regions `heteropolymer_errors.sh`

From reference genomes, computes distribution of heteropolymers depending on the category, i.e. the repeated 2-bases.
Computes a plot (stored in Output/Heteropolymer_errors/heteropolymer_genomic_distribution.png) and raw results containing raw tab-separated table and a LaTeX template of these data (stored in Analysis/$

From explicit alignment files, computes:
 - abundance and error rates depending on heteropolymer category (heteropolymer_abundance_error_rate*)
 - differences between expected heteropolymer length (in reference genome) and sequenced one (difference_expected_sequenced_heteropolymer_length*). Results are displayed as boxplots for which each box merges results for all analysed specie


### (4) Errors associated to trinacleotides regions `trinucleotides.sh`

From reference genomes, computes distribution of trinucleotides depending on the category, i.e. the repeated 3-bases.
Computes occurrences and sequencing accuracy. Results are stored in Output/Trinucleotides and Analysis/Trinucleotides.

---

- [ ] Coming soon: generic version (for the moment low/high GC bact + human) only for low/high GC (remove all "bacteria"/"human" labels)
- [ ] Coming soon: generate a PDF of all computed results, as a summary

---

# Contact

Clara Delahaye: <clara.delahaye@irisa.fr>
