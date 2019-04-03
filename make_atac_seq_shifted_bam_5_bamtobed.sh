#!/bin/bash
#  Aligning fastq files to the genome and creating shifted_sorted.bam file -- preprocessing for ATAC-seq
#  ahead of peak calling
#  
#  STEP 5 of 6
#
#  Name   :  make_atac_seq_shifted_bam_5_bamtobed.sh
#
#  Purpose: This routine sorts the shifted sam file and converts it to a bed file.
#
#  Author:
#            Asli Uyar, PhD
#            Shane Sanders, PhD
#            Anne Deslattes Mays, PhD
#
#  Date:     2016 October 10
#
#  Call:     make_atac_seq_shifted_bam_5_bamtobed.sh
#
#  Assumptions:  
#
#            1.  code has been checked out - fastqc has been run - working directory created
#            2.  we are in the same directory as the scripts.
#            3.  trimmomatic has been run
#            4.  bwa has been run
#            5.  shifted_bam has been called to make the sam file correct for ATAC_seq
#           
#
#-------------------------------------------------------------------------------------
#!/bin/bash

#dataDIR=$1 #ARGV, contains folder with FASTQ files, right now needs full path

scriptDIR=$(pwd)
workingDIR=$scriptDIR/working
inputDIR=$workingDIR/trimmomatic/bwa
outputDIR=$workingDIR/trimmomatic/bwa
#
# now bamtobed
#
rm $workingDIR/shiftedbamfilelist.txt
rm $workingDIR/postPost.qsub

ls -1 $inputDIR/*shifted.bam > $workingDIR/shiftedbamfilelist.txt
FILENUMBER=$(wc -l $workingDIR/shiftedbamfilelist.txt | cut -d' ' -f1)

echo \#!/bin/bash > $workingDIR/postPostBWA.qsub
echo \#PBS -l nodes=1:ppn=16 >> $workingDIR/postPostBWA.qsub
echo \#PBS -l walltime=12:00:00 >> $workingDIR/postPostBWA.qsub
echo \#PBS -N bamtobed  >> $workingDIR/postPostBWA.qsub
echo \#PBS -t 1-$FILENUMBER >> $workingDIR/postPostBWA.qsub
echo module load python >> $workingDIR/postPostBWA.qsub
echo module load R >> $workingDIR/postPostBWA.qsub
echo module load perl/5.16.3 >> $workingDIR/postPostBWA.qsub
echo module load samtools/1.2 >> $workingDIR/postPostBWA.qsub
echo module load bedtools >> $workingDIR/postPostBWA.qsub
echo FILE=\$\(head -n \$PBS_ARRAYID $workingDIR/shiftedbamfilelist.txt \| tail -1\) >> $workingDIR/postPostBWA.qsub
echo FILENAME=\$\(basename \"\${FILE}\" \| sed \'s/\.bam/_sorted.bam/g\'\)  >> $workingDIR/postPostBWA.qsub
echo FILENAME1=\$\(basename \"\${FILE}\" \| sed \'s/\.bam/_sorted/g\'\)  >> $workingDIR/postPostBWA.qsub
echo FILENAME2=\$\(basename \"\${FILE}\" \| sed \'s/\.bam/_sorted.bed/g\'\) >> $workingDIR/postPostBWA.qsub
echo samtools sort \$FILE $outputDIR/\$FILENAME1 >> $workingDIR/postPostBWA.qsub
echo samtools index $outputDIR/\$FILENAME >> $workingDIR/postPostBWA.qsub
echo bedtools bamtobed -i $outputDIR/\$FILENAME \> $outputDIR/\$FILENAME2 >> $workingDIR/postPostBWA.qsub

qsub -V $workingDIR/postPostBWA.qsub

echo "use qstat -u <username> to check the status of your job" 
echo "done with STEP 5 when qsub is complete go to STEP 6 make_atac_seq_shifted_bam_6_cleanup"

