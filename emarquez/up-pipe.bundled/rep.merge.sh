#!/bin/bash

#PBS -w /data/emarquez/
#PBS -l nodes=1:ppn=1,mem=2gb,walltime=24:00:00 
#PBS -m e 
#PBS -M Eladio.Marquez@jax.org 
#PBS -e /home/emarquez/scripts/fperrfile.txt
#PBS -o /dev/null

cd $PBS_O_WORKDIR

module load samtools/1.2

# This scripts simply finds technical replicates for a sample in helix, whose specific locations will vary, merge them and save them in the dedicated directory _bamfiles_ in my data space

outdir="/data/emarquez/bamfiles/HYHO";
indir="/data/auyar/ATAC-Seq_Data/Palucka/HY";

cd ${indir}

samtools merge ${outdir}/HY120_PBMC_50_sorted_merged.bam PBMC/HY120_PBMC_50/HY120_PBMC_50_1/trimmomatic/bwa/HY120-PBMC-50-1_S3_LaneALL_hg19_sorted_rmdup_shifted_sorted.bam PBMC/HY120_PBMC_50/HY120_PBMC_50_2/trimmomatic/bwa/HY120-PBMC-50-2_S4_LaneALL_hg19_sorted_rmdup_shifted_sorted.bam
