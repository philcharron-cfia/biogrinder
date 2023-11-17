# Biogrinder

A versatile omics shotgun and amplicon sequencing read simulator

## Note

I created this repo based off the original biogrinder source code from
https://sourceforge.net/projects/biogrinder/files/biogrinder/.

The tool was originally written in Perl. I have migrated the tool to Python and 
have added some features such as allowing multiple primer sets to be used when
running the amplicon sequencing simulator.

Most of the information found in this file have been adapted from 
https://github.com/zyxue/biogrinder/blob/master/README.md.

## Description

Biogrinder is a versatile program to create random shotgun and amplicon
sequence libraries based on DNA, RNA or proteic reference sequences
provided in a FASTA file.

Biogrinder can produce genomic, metagenomic, transcriptomic,
metatranscriptomic, proteomic, metaproteomic shotgun and amplicon
datasets from current sequencing technologies such as Sanger, 454,
Illumina. These simulated datasets can be used to test the accuracy of
bioinformatic tools under specific hypothesis, e.g. with or without
sequencing errors, or with low or high community diversity. Biogrinder may
also be used to help decide between alternative sequencing methods for a
sequence-based project, e.g. should the library be paired-end or not,
how many reads should be sequenced.

Biogrinder features include:

* shotgun or amplicon read libraries
* omics support to generate genomic, transcriptomic, proteomic, metagenomic,
  metatranscriptomic or metaproteomic datasets
* arbitrary read length distribution and number of reads
* simulation of PCR and sequencing errors (chimeras, point mutations,
  homopolymers)
* support for paired-end (mate pair) datasets
* specific rank-abundance settings or manually given abundance for each genome,
  gene or protein
* creation of datasets with a given richness (alpha diversity)
* independent datasets can share a variable number of genomes (beta diversity)
* modeling of the bias created by varying genome lengths or gene copy number
* profile mechanism to store preferred options
* available to biologists or power users through multiple interfaces: GUI, CLI
  and API

Briefly, given a FASTA file containing reference sequence (genomes,
genes, transcripts or proteins), Biogrinder performs the following steps:

1. Read the reference sequences, and for amplicon datasets, extracts full-length
   reference PCR amplicons using the provided degenerate PCR primers.

1. Determine the community structure based on the provided alpha diversity
   (number of reference sequences in the library), beta diversity (number of
   reference sequences in common between several independent libraries) and
   specified rank- abundance model.

1. Take shotgun reads from the reference sequences or amplicon reads from the
   full- length reference PCR amplicons. The reads may be paired-end reads when
   an insert size distribution is specified. The length of the reads depends on
   the provided read length distribution and their abundance depends on the
   relative abundance in the community structure. Genome length may also biases
   the number of reads to take for shotgun datasets at this step. Similarly, for
   amplicon datasets, the number of copies of the target gene in the reference
   genomes may bias the number of reads to take.

1. Alter reads by inserting sequencing errors (indels, substitutions and
   homopolymer errors) following a position-specific model to simulate reads
   created by current sequencing technologies (Sanger, 454, Illumina). Write the
   reads and their quality scores in FASTA, QUAL and FASTQ files.

## Author

Phil Charron \<<phil.charron@inspection.gc.ca>\>

Adapted from Florent Angly \<<florent.angly@gmail.com>\>

## Installation

Dependencies are found in the [requirements.txt](requirements.txt) file. They 
are installed automatically.

To install Biogrinder globally on your system, run the following commands
in a terminal:

##### On Linux, Unix, MacOS:

```
git clone https://github.com/philcharron-cfia/biogrinder.git
cd biogrinder
python -m pip install --upgrade pip
pip install .
# Test the installation
python tests/run_tests.py
```

## Running Biogrinder

You can run Biogrinder using then command-line interface (CLI). see 
``biogrinder --help```

## Examples

An amplicon library using typical Illumina error rates

```
biogrinder --reference_file 16Sgenes.fna --forward_reverse 16Sprimers.fna --profile_file profile_illumina_amplicon.txt
```

A shotgun DNA library with a coverage of 0.1X

```
biogrinder --reference_file genomes.fna --coverage_fold 0.1
```

Same thing but save the result files in a specific folder and with a specific
name

```
biogrinder --reference_file genomes.fna --coverage_fold 0.1 --base_name my_name --output_dir my_dir
```

A DNA shotgun library with 1000 reads

```
biogrinder --reference_file genomes.fna --total_reads 1000
```

A DNA shotgun library where species are distributed according to a power law

```
biogrinder --reference_file genomes.fna --abundance_model powerlaw 0.1
```

A DNA shotgun library with 123 genomes taken random from the given genomes

```
biogrinder --reference_file genomes.fna --diversity 123
```

Two DNA shotgun libraries that have 50% of the species in common

```
biogrinder --reference_file genomes.fna --num_libraries 2 --shared_perc 50
```

Two DNA shotgun library with no species in common and distributed according to a
exponential rank-abundance model. Note that because the parameter value for the
exponential model is omitted, each library uses a different randomly chosen
value:

```
biogrinder --reference_file genomes.fna --num_libraries 2 --abundance_model exponential
```

A DNA shotgun library where species relative abundances are manually specified

```
biogrinder --reference_file genomes.fna --abundance_file my_abundances.txt
```

A DNA shotgun library with Sanger reads

```
biogrinder --reference_file genomes.fna --read_dist 800 --mutation_dist linear 1 2 --mutation_ratio 80 20
```

A DNA shotgun library with first-generation 454 reads

```
biogrinder --reference_file genomes.fna --read_dist 100 normal 10 --homopolymer_dist balzer
```

A paired-end DNA shotgun library, where the insert size is normally distributed
around 2.5 kbp and has 0.2 kbp standard deviation

```
biogrinder --reference_file genomes.fna --insert_dist 2500 normal 200
```

A transcriptomic dataset

```
biogrinder --reference_file transcripts.fna
```

A unidirectional transcriptomic dataset

```
biogrinder --reference_file transcripts.fna --unidirectional 1
```

A proteomic dataset

```
biogrinder --reference_file proteins.faa --unidirectional 1
```

A 16S rRNA amplicon library
```
biogrinder --reference_file 16Sgenes.fna --forward_reverse 16Sprimers.fna --length_bias 0 --unidirectional 1
```

**Note:** the use of `-length_bias 0` because reference sequence length should
not affect the relative abundance of amplicons.

The same amplicon library with 20% of chimeric reads (90% bimera, 10% trimera)

```
biogrinder --reference_file 16Sgenes.fna --forward_reverse 16Sprimers.fna --length_bias 0 --unidirectional 1 --chimera_perc 20 --chimera_dist 90 10
```

Three 16S rRNA amplicon libraries with specified MIDs and no reference sequences
in common

```
biogrinder --reference_file 16Sgenes.fna --forward_reverse 16Sprimers.fna --length_bias 0 --unidirectional 1 --num_libraries 3 --multiplex_ids MIDs.fna
```

Reading reference sequences from the standard input, which allows you to
decompress FASTA files on the fly:

```
zcat microbial_db.fna.gz | biogrinder --reference_file - --total_reads 100
```

## CLI Output

For each shotgun or amplicon read library requested, the following files are
generated:

1. A rank-abundance file, tab-delimited, that shows the relative abundance of
   the different reference sequences.
1. A library summary file that shows details about the created library.
1. If amplicon method is chosen rathen than shotguen, an amplicon summary file 
   with details about amplicons that were identified.
1. A file containing the read sequences in FASTA format. The read headers
   contain information necessary to track from which reference sequence each
   read was taken and what errors it contains. This file is not generated if
   `fastq_output` option was provided.
1. If the `qual_levels` option was specified, a file containing the quality
   scores of the reads (in QUAL format).
1. If the `fastq_output` option was provided, a file containing the read
   sequences in FASTQ format.

## Features in Future Versions
- Mutation rates and documentation to allow long-read simulation
- Allow mismatches in primers when using amplicon sequencing

## Copyright

Copyright 2023 Phil Charron (phil.charron@inspection.gc.ca)

Biogrinder is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License (GPL) as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version. Biogrinder is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with
Biogrinder. If not, see http://www.gnu.org/licenses/.


## Bugs

If you have any issues installing or running VariantDetective, or would like a new feature added to the tool, please open an issue here on GitHub.