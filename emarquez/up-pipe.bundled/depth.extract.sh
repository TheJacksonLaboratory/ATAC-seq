#!/bin/bash

#PBS -w /data/emarquez/peaks
#PBS -l nodes=1:ppn=1,mem=2gb,walltime=48:00:00 
#PBS -m abe 
#PBS -M Eladio.Marquez@jax.org 
#PBS -e /home/emarquez/scripts/fperrfile.txt
#PBS -o /dev/null

module load samtools/1.2

# Compiles depth of all samples into a single file

> /data/emarquez/peaks/bamcounts/pooledsamp_depth.txt.tmp;
for d in Dir1 Dir2 ...;
 do
  for f in $(find /data/emarquez/bamfiles/${d}/ -maxdepth 1 -name '*.bam');
   do
    export name=$f;
    samtools view -c ${f} | perl -lane 'print "$ENV{name},",$F[0]' >> /data/emarquez/peaks/bamcounts/pooledsamp_depth.txt.tmp;
 done;
done;

sort -k1 /data/emarquez/peaks/bamcounts/pooledsamp_depth.txt.tmp | uniq > /data/emarquez/peaks/bamcounts/blah_pooledsamp_depth.csv;
rm /data/emarquez/peaks/bamcounts/pooledsamp_depth.txt.tmp;