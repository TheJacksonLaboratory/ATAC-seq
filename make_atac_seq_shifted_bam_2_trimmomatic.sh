#!/bin/bash
#  Aligning fastq files to the genome and creating shifted_sorted.bam file -- preprocessing for ATAC-seq
#  ahead of peak calling
#  
#  STEP 2 of 6
#
#  Name   :  make_atac_seq_shifted_bam_2_trimmomatic.sh
#
#  Purpose:  fastqc has been run now we run trimmomatic to identified the trimmend and untrimmed fastq files
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
#
#-------------------------------------------------------------------------------------
dataDIR=$1 #ARGV, contains folder with FASTQ files, right now needs trailing /

scriptDIR=$(pwd)
outputDIR=$scriptDIR/working

## Trimmomatic & Adapter Trimming Pipeline ###
rm $outputDIR/trimmomatic_R1_filelist.txt
rm $outputDIR/trimmomatic.qsub

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
echo FILE1paired=\$\(basename \"\${FILE}\" \| sed \'s/R1_001\.fastq.gz/R1_001\.trim\.fastq\.gz/g\'\) >> $outputDIR/trimmomatic.qsub
echo FILE2paired=\$\(basename \"\${FILE2}\" \| sed \'s/R2_001\.fastq.gz/R2_001\.trim\.fastq\.gz/g\'\) >> $outputDIR/trimmomatic.qsub
echo FILE1unpaired=\$\(basename \"\${FILE}\" \| sed \'s/R1_001\.fastq.gz/R1_001\.trimU\.fastq\.gz/g\'\) >> $outputDIR/trimmomatic.qsub
echo FILE2unpaired=\$\(basename \"\${FILE2}\" \| sed \'s/R2_001\.fastq.gz/R2_001\.trimU\.fastq\.gz/g\'\) >> $outputDIR/trimmomatic.qsub
echo java -jar /opt/compsci/Trimmomatic/0.33/trimmomatic-0.33.jar PE -threads 2 \$FILE $dataDIR\$FILE2 $outputDIR/trimmomatic/\$FILE1paired $outputDIR/trimmomatic/\$FILE1unpaired $outputDIR/trimmomatic/\$FILE2paired $outputDIR/trimmomatic/\$FILE2unpaired TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 >> $outputDIR/trimmomatic.qsub
echo python $scriptDIR/auyar/pyadapter_trim.py -a $outputDIR/trimmomatic/\$FILE1paired -b $outputDIR/trimmomatic/\$FILE2paired >> $outputDIR/trimmomatic.qsub

######
echo about to call $outputDIR/trimmomatic.qsub
qsub -V $outputDIR/trimmomatic.qsub

echo "made the call to $outputDIR/trimmomatic.qsub -- wait until that job is complete before going to STEP 3"
echo "use qstat -u <username> to check the status of your job" 
echo "done with STEP 2 when qsub is complete go to STEP 3 make_atac_seq_shifted_bam_3_bwa"
