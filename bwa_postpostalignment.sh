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
rm $workingDIR/postPostBWA.qsub

ls -1 $inputDIR/*shifted.bam > $workingDIR/shiftedbamfilelist.txt
FILENUMBER=$(wc -l $workingDIR/shiftedbamfilelist.txt | cut -d' ' -f1)

echo \#!/bin/bash > $workingDIR/postPostBWA.qsub
echo \#PBS -l nodes=1:ppn=16 >> $workingDIR/postPostBWA.qsub
echo \#PBS -l walltime=12:00:00 >> $workingDIR/postPostBWA.qsub
echo \#PBS -N bwa  >> $workingDIR/postPostBWA.qsub
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

