#!/bin/bash

#PBS -l nodes=1:ppn=12,walltime=16:00:00
#PBS -q batch
  

cd $PBS_O_WORKDIR

module load python
module load R 
module load perl/5.10.1
module load samtools/0.1.19
module load bedtools

#gunzip ../*
#mkdir ../fastQC
#mkdir ../trimmomatic

#cd ../trimmomatic
#mkdir bwa
#cd bwa
#mkdir macs2
#mkdir homer
#cd ../..

#cd scripts
#bash fastQC.sh ../
#bash trimmomatic.sh ../ 
#mv *.trim.fastq ../trimmomatic/

bash bwa_trimmomatic.sh ../
bash macs2_bwa_trimmomatic.sh ../
bash homer_bwa_trimmomatic.sh ../
bash nReads.sh ../





