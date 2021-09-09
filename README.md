# QASeq
CNV analysis code for QASeq manuscript

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Part 1: obtaining UMI family counts from NGS fastq files

Input: R1 151nt + R2 151nt paired-end fastq file; example: "lib1_R1.fastq" and "lib1_R2.fastq"

Output: Spreadsheet containing UMI family counts of each amplicon: "allcountsDynamic.csv"

Code:

shell_commands.sh: Main MacOS shell script running all the analysis code

adapter_trim_PE_downsamp_lrw_20210831.py: Python code for Illumina adapter trimming, can perform downsampling of reads

sort_index.py: Python code for sorting SAM file

UMIcounter.py: Python code for UMI grouping

analysis_umisorted_v2_Dynamic.m: Matlab code for UMI counting

Other files:

all2seq_new.fasta, all179seq.fasta, all226seq.fasta: Amplicon sequences of different panels

Instructions for use:

In terminal, navigate to the current folder containing codes and the fastq file to analyze, then type in "sh shell_commands.sh".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Part 2: CNV calling based on sequential Mann-Whitney U test

Input: UMI family counts and plex information, example "example_UMI_count_input.xlsx"

Output: Utest results and figures in Matlab

Code:

example_Utest_plot.m: Matlab code for sequential Mann-Whitney U test

UtestQASeq_allPvalv2.m: U test core function

Other files:

example_UMI_count_input.xlsx: example input

All_StdSamples.xlsx: molecule counts of standard (healthy) samples for conversion yield calculation

Instructions for use:

Place the required filed of Part 2 in the same folder; open Matlab and run "example_Utest_plot.m".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Requirement:

The software has been tested on MacOS (10.13.3 or 10.15.7) with the following installed:

Python 2.7.16 or Python 2.7.18

Matlab R2019a or Matlab R2020b

Bowtie 2 version 2.3.5.1

pysam
