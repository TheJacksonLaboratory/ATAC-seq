#!/bin/bash
#  Rename file names
#
#  Name   :  make_atac_seq_rename.sh
#
#  Purpose:  Use more readable names for samples. 
#
#  Author:
#            Duygu Ucar, PhD
#
#  Date:     2017 November 15
#
#  Call:     sh make_atac_seq_rename.sh <tab delimited file with SRA id sampleID> <folder containing macs output>
#
#  Assumptions:  
#
#            1.  ATAC-seq processing is done
#            2.  Peaks are called
#
#-------------------------------------------------------------------------------------
filename=$1 # ARGV, file that lists the new names
peakDIR=$2 #ARGV, contains folder with peak call results right now needs trailing /

scriptDIR=$(pwd)

#while read -r line
while read sra label
do
    echo "$sra:$label:"
    FILENAME=$sra"_1.trim.trim.fastq.hg19_sorted_rmdup_shifted_sorted"
    rename $FILENAME $label $peakDIR/*
    echo $FILENAME
done < "$filename"

