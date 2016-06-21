#!/bin/bash

dataDIR=$1

outputDIR=/home/ssander/Desktop/ATAC-seq/working
mkdir $outputDIR

### FastQC Pipeline
mkdir $outputDIR/fastQC

find $dataDIR -name *.fastq.gz > $outputDIR/fastqc_filelist.txt
FILENUMBER=$(wc -l $outputDIR/fastqc_filelist.txt | cut -d' ' -f1)

echo $FILENUMBER

echo \#!/bin/bash >> $outputDIR/fastqc.qsub
echo \#PBS -l nodes=1:ppn=2 >> $outputDIR/fastqc.qsub
echo \#PBS -l walltime=24:00:00 >> $outputDIR/fastqc.qsub
echo \#PBS -N fastqc  >> $outputDIR/fastqc.qsub
echo \#PBS -t 1-$FILENUMBER >> $outputDIR/fastqc.qsub
echo FILE=\$\(head -n \$PBS_ARRAYID $outputDIR/fastqc_filelist.txt \| tail -1\) >> $outputDIR/fastqc.qsub
echo /opt/compsci/FastQC/0.11.3/fastqc -t 2 --noextract \$FILE -o $outputDIR/fastQC >> $outputDIR/fastqc.qsub 

## Trimmomatic 
mkdir $outputDIR/trimmomatic
find $dataDIR -name *R1_001.fastq.gz > $outputDIR/trimmomatic_R1_filelist.txt
FILENUMBER=$(wc -l $outputDIR/trimmomatic_R1_filelist.txt | cut -d' ' -f1)

echo $FILENUMBER

echo \#!/bin/bash >> $outputDIR/trimmomatic.qsub
echo \#PBS -l nodes=1:ppn=2 >> $outputDIR/trimmomatic.qsub
echo \#PBS -l walltime=24:00:00 >> $outputDIR/trimmomatic.qsub
echo \#PBS -N fastqc  >> $outputDIR/trimmomatic.qsub
echo \#PBS -t 1-$FILENUMBER >> $outputDIR/trimmomatic.qsub
echo FILE=\$\(head -n \$PBS_ARRAYID $outputDIR/fastqc_filelist.txt \| tail -1\) >> $outputDIR/trimmomatic.qsub
echo FILE2=\$\(basename "\${FILE}" \| sed \'s/R1_001\.fastq.gz/R2_001\.fastq.gz/g\'\) >> $outputDIR/trimmomatic.qsub
# need to finish adding auyuar/trimmomatic.sh here #


#qsub $outputDIR/fastqc.qsub
#qsub $outputDIR/trimmomatic.sh
