#!/bin/bash

dataDIR=$1 #ARGV, contains folder with FASTQ files, right now needs trailing /

scriptDIR=$(pwd)
outputDIR=$scriptDIR/working
mkdir $outputDIR

### FastQC Pipeline ###
mkdir $outputDIR/fastQC

find $dataDIR -name *.fastq.gz > $outputDIR/fastqc_filelist.txt
FILENUMBER=$(wc -l $outputDIR/fastqc_filelist.txt | cut -d' ' -f1)

#echo $FILENUMBER

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

#echo $FILENUMBER

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
echo python $scriptDIR/auyar/pyadapter_trim.py -a $outputDIR/trimmomatic/\$FILE1paired -b $outputDIR/trimmomatic/\$FILE2paired >> $outputDIR/trimmomatic.qsub

######

### ATAC-seq Portion ###
###
### This section is still a work in progress...
###
###
# Next steps...
# auyar/run_ATAC_pipeline.sh, which calls:
# * auyar/bwa_trimmomatic.sh
# * auyar/macs2_bwa_trimmomatic.sh
# * auyar/homer_bwa_trimmomatic.sh
# * auyar/nReads.sh

cd $outputDIR/trimmomatic
mkdir bwa
cd bwa
mkdir macs2
mkdir homer
cd ../..
currentDIR=$(pwd)
#echo $currentDIR

echo \#!/bin/bash >> $outputDIR/ATAC-seq.qsub
echo \#PBS -l nodes=1:ppn=16 >> $outputDIR/ATAC-seq.qsub
echo \#PBS -l walltime=16:00:00 >> $outputDIR/ATAC-seq.qsub
echo \#PBS -N ATAC-seq  >> $outputDIR/ATAC-seq.qsub
echo module load python >> $outputDIR/ATAC-seq.qsub
echo module load R >> $outputDIR/ATAC-seq.qsub
echo module load perl/5.10.1 >> $outputDIR/ATAC-seq.qsub
echo module load samtools/0.1.19 >> $outputDIR/ATAC-seq.qsub
echo module load bedtools >> $outputDIR/ATAC-seq.qsub
echo gunzip $outputDIR/trimmomatic/*.gz >> $outputDIR/ATAC-seq.qsub
#echo bash $scriptDIR/auyar/bwa_trimmomatic.sh ../ >> $outputDIR/ATAC-seq.qsub

ls $outputDIR/trimmomatic
for f in $(find $outputDIR/trimmomatic/ -name *R1_001.trim.fastq.gz)
do
	fileName2=$(basename "${f}" | sed 's/R1_001\.trim\.fastq/R2_001\.trim\.fastq/g')
	fileNameSAM=$(basename "${f}" | sed 's/R1_001\.trim\.fastq/hg19\.sam/g')
	fileNameMetrics=$(basename "${f}" | sed 's/R1_001\.trim\.fastq/metrics\.txt/g')
	echo /opt/compsci/bwa/0.7.12/bin/bwa mem -M /data/shared/genomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa $f $outputDIR/trimmomatic/$fileName2 > $outputDIR/trimmomatic/bwa/$fileNameSAM
done


echo bash $scriptDIR/auyar/macs2_bwa_trimmomatic.sh ../ >> $outputDIR/ATAC-seq.qsub
echo bash $scriptDIR/auyar/homer_bwa_trimmomatic.sh ../ >> $outputDIR/ATAC-seq.qsub
echo bash $scriptDIR/auyar/nReads.sh ../ >> $outputDIR/ATAC-seq.qsub
echo gzip $outputDIR/trimmomatic/*.fastq >> $outputDIR/ATAC-seq.qsub
######

# something buggy is going on here...do job dependencies bug out when given job arrays?
# commenting the automatic queue submission section out for debugging at this point
# qsub $outputDIR/fastqc.qsub
# JOBHOLD="$(qsub $outputDIR/trimmomatic.qsub)"
# qsub -Wdepend=afterok:$JOBHOLD $outputDIR/ATAC-seq.qsub
