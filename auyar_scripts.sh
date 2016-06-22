#!/bin/bash

scriptDIR=$(pwd)

dataDIR=$1

outputDIR=/home/ssander/Desktop/ATAC-seq/working
mkdir $outputDIR

### FastQC Pipeline ###
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
######

## Trimmomatic & Adapter Trimming Pipeline ###
mkdir $outputDIR/trimmomatic
find $dataDIR -name *R1_001.fastq.gz > $outputDIR/trimmomatic_R1_filelist.txt
FILENUMBER=$(wc -l $outputDIR/trimmomatic_R1_filelist.txt | cut -d' ' -f1)

echo $FILENUMBER

echo \#!/bin/bash >> $outputDIR/trimmomatic.qsub
echo \#PBS -l nodes=1:ppn=2 >> $outputDIR/trimmomatic.qsub
echo \#PBS -l walltime=24:00:00 >> $outputDIR/trimmomatic.qsub
echo \#PBS -N trimmomatic  >> $outputDIR/trimmomatic.qsub
echo \#PBS -t 1-$FILENUMBER >> $outputDIR/trimmomatic.qsub
echo module load python >> $outputDIR/trimmomatic.qsub
echo FILE=\$\(head -n \$PBS_ARRAYID $outputDIR/trimmomatic_R1_filelist.txt \| tail -1\) >> $outputDIR/trimmomatic.qsub
echo FILE2=\$\(basename "\${FILE}" \| sed \'s/R1_001\.fastq.gz/R2_001\.fastq.gz/g\'\) >> $outputDIR/trimmomatic.qsub
echo FILE1paired=\$\(basename \"\${FILE}\" \| sed \'s/R1_001\.fastq.gz/R1_001\.fastq_filtered.gz/g\'\) >> $outputDIR/trimmomatic.qsub
echo FILE2paired=\$\(basename \"\${FILE2}\" \| sed \'s/R2_001\.fastq.gz/R2_001\.fastq_filtered.gz/g\'\) >> $outputDIR/trimmomatic.qsub
echo FILE1unpaired=\$\(basename \"\${FILE}\" \| sed \'s/R1_001\.fastq.gz/R1_001\.trimU\.fastq.gz/g\'\) >> $outputDIR/trimmomatic.qsub
echo FILE2unpaired=\$\(basename \"\${FILE2}\" \| sed \'s/R2_001\.fastq.gz/R2_001\.trimU\.fastq.gz/g\'\) >> $outputDIR/trimmomatic.qsub
echo java -jar /opt/compsci/Trimmomatic/0.33/trimmomatic-0.33.jar PE -threads 2 \$FILE $dataDIR\$FILE2 $outputDIR/trimmomatic/\$FILE1paired $outputDIR/trimmomatic/\$FILE1unpaired $outputDIR/trimmomatic/\$FILE2paired $outputDIR/trimmomatic/\$FILE2unpaired TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 >> $outputDIR/trimmomatic.qsub

# Problems with the pyadapter_trim.py script and dependencies below:

echo python $scriptDIR/auyar/pyadapter_trim.py -a $outputDIR/trimmomatic/\$FILE1paired -b $outputDIR/trimmomatic/\$FILE2paired >> $outputDIR/trimmomatic.qsub

######

# Next steps...
# auyar/run_ATAC_pipeline.sh, which calls:
# * auyar/bwa_trimmomatic.sh
# * auyar/macs2_bwa_trimmomatic.sh
# * auyar/homer_bwa_trimmomatic.sh
# * auyar/nReads.sh
######

#qsub $outputDIR/fastqc.qsub
#qsub $outputDIR/trimmomatic.sh
