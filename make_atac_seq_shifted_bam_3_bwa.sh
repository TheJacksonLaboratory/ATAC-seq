#!/bin/bash
#  Aligning fastq files to the genome and creating shifted_sorted.bam file -- preprocessing for ATAC-seq
#  ahead of peak calling
#  
#  STEP 3 of 6
#
#  Name   :  make_atac_seq_shifted_bam_3_bwa.sh
#
#  Purpose:  bwa is the alignment of the trimmed fastq files to the genome
#             
#  Author:
#            Asli Uyar, PhD
#            Shane Sanders, PhD
#            Anne Deslattes Mays, PhD
#
#  Date:     2016 October 10
#
#  Call:     make_atac_seq_shifted_bam_1_run_fastqc.sh <ATAC-Seq Banchereau-Lab/GT-delivery/ATAC-seq directory (with trailing /> ]
#
#  Assumptions:  
#
#            1.  code has been checked out - fastqc has been run - working directory created
#            2.  we are in the same directory as the scripts.
#            3.  trimmomatic has been run
#
#-------------------------------------------------------------------------------------
#dataDIR=$1 #ARGV, contains folder with FASTQ files, right now needs trailing /

scriptDIR=$(pwd)
workingDIR=$scriptDIR/working

## Bwa Alignment Pipeline ##
mkdir $workingDIR/trimmomatic/bwa

rm $workingDIR/filelist.txt
ls -1 $workingDIR/*R1_001.trim.fastq.gz > $workingDIR/filelist.txt
FILENUMBER=$(wc -l $workingDIR/filelist.txt | cut -d' ' -f1)

#echo $FILENUMBER
rm $workingDIR/bwa.qsub

echo \#!/bin/bash >> $workingDIR/bwa.qsub
echo \#PBS -l nodes=1:ppn=16 >> $workingDIR/bwa.qsub
echo \#PBS -l walltime=48:00:00 >> $workingDIR/bwa.qsub
echo \#PBS -N bwa  >> $workingDIR/bwa.qsub
echo \#PBS -t 1-$FILENUMBER%100 >> $workingDIR/bwa.qsub
echo module load python >> $workingDIR/bwa.qsub
echo module load R >> $workingDIR/bwa.qsub
echo module load perl/5.10.1 >> $workingDIR/bwa.qsub
echo module load samtools/0.1.19 >> $workingDIR/bwa.qsub
echo module load bedtools >> $workingDIR/bwa.qsub
echo FILE=\$\(head -n \$PBS_ARRAYID $workingDIR/filelist.txt \| tail -1\) >> $workingDIR/bwa.qsub
echo FILE2=\$\(basename "\${FILE}"\| sed \'s/R1_001\.trim\.fastq\.gz/R2_001\.trim\.fastq\.gz/g\'\) >> $workingDIR/bwa.qsub
echo FILESAM=\$\(basename "\${FILE}"\).hg19.sam >> $workingDIR/bwa.qsub
echo /opt/compsci/bwa/0.7.12/bin/bwa mem -M /data/shared/genomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa \$FILE $workingDIR/\$FILE2 \> $workingDIR/bwa/\$FILESAM >> $workingDIR/bwa.qsub

######
qsub -V $workingDIR/bwa.qsub

echo "made the call to $workingDIR/bwa.qsub -- wait until that job is complete before going to STEP 3"
echo "use qstat -u <username> to check the status of your job" 
echo "done with STEP 3 when qsub is complete go to STEP 4 make_atac_seq_shifted_bam_4_bwa_bam_shift"

