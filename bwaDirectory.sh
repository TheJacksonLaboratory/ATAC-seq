#!/bin/bash

dataDIR=$1 #ARGV, contains folder with FASTQ files, right now needs full path
rm $dataDIR/filelist.txt
#find $dataDIR -name *R1_001.trim.fastq.gz > $dataDIR/filelist.txt
ls -1 $dataDIR/*R1_001.trim.fastq.gz > $dataDIR/filelist.txt
FILENUMBER=$(wc -l $dataDIR/filelist.txt | cut -d' ' -f1)

#echo $FILENUMBER
rm $dataDIR/bwa.qsub

echo \#!/bin/bash >> $dataDIR/bwa.qsub
echo \#PBS -l nodes=1:ppn=16 >> $dataDIR/bwa.qsub
echo \#PBS -l walltime=48:00:00 >> $dataDIR/bwa.qsub
echo \#PBS -N bwa  >> $dataDIR/bwa.qsub
echo \#PBS -t 1-$FILENUMBER%100 >> $dataDIR/bwa.qsub
echo module load python >> $dataDIR/bwa.qsub
echo module load R >> $dataDIR/bwa.qsub
echo module load perl/5.10.1 >> $dataDIR/bwa.qsub
echo module load samtools/0.1.19 >> $dataDIR/bwa.qsub
echo module load bedtools >> $dataDIR/bwa.qsub
echo FILE=\$\(head -n \$PBS_ARRAYID $dataDIR/filelist.txt \| tail -1\) >> $dataDIR/bwa.qsub
echo FILE2=\$\(basename "\${FILE}"\| sed \'s/R1_001\.trim\.fastq\.gz/R2_001\.trim\.fastq\.gz/g\'\) >> $dataDIR/bwa.qsub
echo FILESAM=\$\(basename "\${FILE}"\).hg19.sam >> $dataDIR/bwa.qsub
echo /opt/compsci/bwa/0.7.12/bin/bwa mem -M /data/shared/genomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa \$FILE $dataDIR/\$FILE2 \> $dataDIR/bwa/\$FILESAM >> $dataDIR/bwa.qsub

######








#
