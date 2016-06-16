#!/bin/bash

#PBS -w /data/emarquez/peaks
#PBS -l nodes=1:ppn=1,mem=8gb,walltime=48:00:00 
#PBS -m abe 
#PBS -M Eladio.Marquez@jax.org 
#PBS -e /home/emarquez/scripts/fperrfile.txt
#PBS -o /dev/null

module load bedtools/2.17.0

# This script retrieves peak signals and scores from peak files at consensus peaks. Also includes peak width and overlapping width for each peak.
# Input is first three columns (chr/start/end) of consensus and each BED file originally used to build the consensus
# Output is a peakStats file for each original sample

cd /data/emarquez/peaks/
peakset=$(find Dir1/*.broadPeak) # Dir1 is the given name of a group of samples

for f in $peakset
 do
 	cut -f1-3 bamcounts/blah_consensus_peaks_rawcounts.txt \
 	| perl -ne 'print if $.>1' | intersectBed -wao -a - -b $f | cut -f1-3,5-8,10-13 \
 	| perl -lane 'print $F[0],"\t",$F[1],"\t",$F[2],"\t",$F[4]-$F[3],"\t",$F[5],"\t",$F[6],"\t",$F[7],"\t",$F[8],"\t",$F[9],"\t",$F[10]' > consensus/Dir1/$(basename $f .broadPeak).peakStats
done;