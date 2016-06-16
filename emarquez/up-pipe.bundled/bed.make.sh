#!/bin/bash

#PBS -w /data/emarquez/peaks
#PBS -l nodes=1:ppn=1,mem=2gb,walltime=148:00:00 
#PBS -m abe 
#PBS -M Eladio.Marquez@jax.org 
#PBS -e /home/emarquez/scripts/fperrfile.txt
#PBS -o /dev/null

module load bedtools/2.17.0

# This script cycles through appropriate BAM files in the bamfiles directory and create BED files for each, which output to dedicated bedfiles directory in my data space

for d in Dir1 Dir2 ...;
 do
  for f in $(find /data/emarquez/bamfiles/${d}/ -maxdepth 1 -name '*.bam');
   do
    name=$(basename "${f}" .bam);
    bamToBed -i /data/emarquez/bamfiles/${d}/${name}.bam > /data/emarquez/bedfiles/${d}/${name}.bed
 done;
done;