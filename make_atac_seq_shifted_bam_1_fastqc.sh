#!/bin/bash
#  Aligning fastq files to the genome and creating shifted_sorted.bam file -- preprocessing for ATAC-seq
#  ahead of peak calling
#  
#  STEP 1 of 6
#
#  Name   :  make_atac_seq_shifted_bam_1_fastqc.sh
#
#  Purpose:  Once it is known that the files are ready from Genome Technologies (GT)
#            These files are copied to /data/Banchereau-Lab/GT-delivery/ATAC-seq maintaining
#            the project information encoded in the directory name
#    
#            Once copied, the user goes to the /data/Banchereau-Lab/ATAC-seq directory
#            makes an appropriately named subdirectory -- and begins the process
#            this may be facilitated by aligning all files in the project at once
#            once they are all delivered for example, and it may be over the entire course
#            of several months -- a single run may be performed and the job can be job arrayed
#            
#            To facilitate this -- symbolic links are made to a representatively named
#            directory structure and these may then be referenced in this first call
#
#  Use Case: Upon completion of all the sequencing that is to be performed
#            A user creates a directory, say hyho_cd8_nk_nkt_depletion_analysis in the
#            /data/Banchereau-Lab/GT-delivery/ATAC-seq directory/hyho_cd8_nk_nkt_depletion_analysis
#            the user makes a symbolic link to all the files that are necessary to run
#            or in this case, the remaining files after part of them have been already run
#            then makes a similarily named file in the /data/Banchereau-Lab/ATAC-seq directoy.
#            git clones the ATAC-seq files and runs this workflow.
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
#
#-------------------------------------------------------------------------------------

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

qsub -V $outputDIR/fastqc.qsub

echo "made the call to $outputDIR/fastqc.qsub -- wait until that job is complete before going to STEP 3"
echo "use qstat -u <username> to check the status of your job" 
echo "done with STEP 1 when qsub is complete go to STEP 2 make_atac_seq_shifted_bam_2_trimmomatic.sh"

