#!/bin/bash

#PBS -w /data/emarquez/peaks
#PBS -l nodes=1:ppn=1,mem=2gb,walltime=48:00:00 
#PBS -m abe 
#PBS -M Eladio.Marquez@jax.org 
#PBS -e /home/emarquez/scripts/fperrfile.txt
#PBS -o /dev/null

module load bedtools/2.17.0

# This script reads peaks from multiple files from directory peaks, combines them via multi-intersect, merge the resulting bookended regions to create non-overlapping consensus peaks.
# The last line subsets the consensus peaks to include only those that have been called in 3 samples or more.

bedtools multiinter -i $(find /data/emarquez/peaks/Dir1/ -maxdepth 1 -name "*.bed") > /data/emarquez/peaks/consensus/blah_sorted_merged_peaks_multiinter.txt;
bedtools merge -i /data/emarquez/peaks/consensus/aging_sorted_merged_peaks_multiinter.txt > /data/emarquez/peaks/consensus/consensus_peaks.tmp;
bedtools annotate -counts -i /data/emarquez/peaks/consensus/consensus_peaks.tmp -files $(find /data/emarquez/peaks/Dir1/ -maxdepth 1 -name "*.bed")	| sortBed | cut -f4- | sed 's/[1-9]/1/g' | sed 's/	/+/g' | sed 's/+$//g' | bc | paste /data/emarquez/peaks/consensus/consensus_peaks.tmp - > /data/emarquez/peaks/consensus/blah_sorted_merged_consensus_peaks.bed;
perl -lane 'print if $F[3]>2' /data/emarquez/peaks/consensus/blah_sorted_merged_consensus_peaks.bed > /data/emarquez/peaks/consensus/blah_sorted_merged_consensus_peaks_minOverlap3.bed;
rm /data/emarquez/peaks/consensus/consensus_peaks.tmp;