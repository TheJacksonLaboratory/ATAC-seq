#!/bin/bash

dataDIR=$1 #ARGV, contains folder with FASTQ files, right now needs full path

scriptDIR=$(pwd)
inputDIR=$dataDIR/trimmomatic/bwa
outputDIR=$dataDIR/trimmomatic/bwa
#
# now bamtobed
#
rm $dataDIR/shiftedbamfilelist.txt
rm $dataDIR/postPost.qsub

ls -1 $inputDIR/*shifted_sorted.bam > $dataDIR/shiftedbamfilelist.txt
FILENUMBER=$(wc -l $dataDIR/shiftedbamfilelist.txt | cut -d' ' -f1)

echo \#!/bin/bash > $dataDIR/postPostBWA.qsub
echo \#PBS -l nodes=1:ppn=16 >> $dataDIR/postPostBWA.qsub
echo \#PBS -l walltime=12:00:00 >> $dataDIR/postPostBWA.qsub
echo \#PBS -N bwa  >> $dataDIR/postPostBWA.qsub
echo \#PBS -t 1-$FILENUMBER >> $dataDIR/postPostBWA.qsub
echo module load python >> $dataDIR/postPostBWA.qsub
echo module load R >> $dataDIR/postPostBWA.qsub
echo module load perl/5.10.1 >> $dataDIR/postPostBWA.qsub
echo module load samtools/0.1.19 >> $dataDIR/postPostBWA.qsub
echo module load bedtools >> $dataDIR/postPostBWA.qsub
echo FILE=\$\(head -n \$PBS_ARRAYID $dataDIR/shiftedbamfilelist.txt \| tail -1\) >> $dataDIR/postPostBWA.qsub
echo FILENAME=\$\(basename \"\${FILE}\" \| sed \'s/\.bam/_sorted.bam/g\'\)  >> $dataDIR/postPostBWA.qsub
echo FILENAME1=\$\(basename \"\${FILE}\" \| sed \'s/\.bam/_sorted/g\'\)  >> $dataDIR/postPostBWA.qsub
echo FILENAME2=\$\(basename \"\${FILE}\" \| sed \'s/\.bam/_sorted.bed/g\'\) >> $dataDIR/postPostBWA.qsub
echo samtools sort \$FILE $outputDIR/\$FILENAME1 >> $dataDIR/postPostBWA.qsub
echo samtools index $outputDIR/\$FILENAME >> $dataDIR/postPostBWA.qsub
echo bedtools bamtobed -i $outputDIR/\$FILENAME \> $outputDIR/\$FILENAME2 >> $dataDIR/postPostBWA.qsub

