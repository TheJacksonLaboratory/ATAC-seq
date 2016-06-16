#!/bin/bash

#PBS -w /data/emarquez/peaks
#PBS -l nodes=1:ppn=1,mem=2gb,walltime=48:00:00 
#PBS -m abe 
#PBS -M Eladio.Marquez@jax.org 
#PBS -e /home/emarquez/scripts/fperrfile.txt
#PBS -o /dev/null

module load samtools/0.1.19

# This scripts created BAI files for merged BAM files and saves them in the same directory

for d in Dir1 Dir2 ...;
 do
 for t in type1 type2 ...; # from naming convention
 do
  for f in $(find /data/emarquez/bamfiles/${d}/ -maxdepth 1 -name "*_${t}_*.bam");
   do
    samtools index ${f};
   done;
 done;
done;