#!/bin/bash
#  Aligning fastq files to the genome and creating shifted_sorted.bam file -- preprocessing for ATAC-seq
#  ahead of peak calling
#  
#  STEP 6 of 6
#
#  Name   :  make_atac_seq_shifted_bam_6_cleanup.sh
#
#  Purpose: This routine removes all the trim and trimU fastq files (intermediate files that may be reproduced)
#           It also removes all sam files, and the unsorted bam file, also all .o and .e files for readibility.
#
#  Author:
#            Asli Uyar, PhD
#            Shane Sanders, PhD
#            Anne Deslattes Mays, PhD
#
#  Date:     2016 October 10
#
#  Call:     make_atac_seq_shifted_bam_6_cleanup.sh <ATAC-Seq Banchereau-Lab/GT-delivery/ATAC-seq directory (with trailing /> ]
#
#  Assumptions:  
#
#            1.  code has been checked out - fastqc has been run - working directory created
#            2.  we are in the same directory as the scripts.
#            3.  trimmomatic has been run
#            4.  bwa has been run
#            5.  shifted_bam has been called to make the sam file correct for ATAC_seq
#            6.  the shifted_bam script and bedtobam has been called so all that is left to do is cleanup!
#           
#
#-------------------------------------------------------------------------------------
#dataDIR=$1 #ARGV, contains folder with FASTQ files, right now needs trailing /

scriptDIR=$(pwd)
workingDIR=$scriptDIR/working
trimmomaticDIR=$workingDIR/trimmomatic
bwaDIR=$trimmmoaticDIR/trimmomatic/bwa

rm $scriptDIR/*.o*
rm $scriptDIR/*.e*

rm $trimmomaticDIR/*.o*
rm $trimmomaticDIR/*.e*

rm $trimmomaticDIR/*.trim.*
rm $trimmomaticDIR/*.trimU.*

rm $bwaDIR/*.sam
rm $bwaDIR/*shifted.bam

echo "done with STEP 6 cleanup is complete all is done with make_atac_seq_shifted_bam"

