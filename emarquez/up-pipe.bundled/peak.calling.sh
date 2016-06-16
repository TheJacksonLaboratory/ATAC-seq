#!/bin/bash

#PBS -w /data/emarquez/tagdir
#PBS -l nodes=1:ppn=1,mem=2gb,walltime=48:00:00 
#PBS -m abe 
#PBS -M Eladio.Marquez@jax.org 
#PBS -e /home/emarquez/scripts/fperrfile.txt
#PBS -o /dev/null

module load Python/2.7.3
module load MACS/2.1.0.20151222

# This script was used to call peaks on merged datasets (BED files) for the aging and other projects. Output is sent to dedicated folder peaks

d=<Dir1>;
g=<agegroup>; # HY, HO, HM, F, VHY
for f in $(find /data/emarquez/bedfiles/${d}/ -maxdepth 1 -name "${g}*.bed");
 do
  echo ${f};
  export thisn=$(basename $f .bed);
  thispeaks="/data/emarquez/peaks/${d}/$(basename ${f} .bed)";
  macs2 callpeak -t $f -f BED -n $thispeaks -g 'hs' --nomodel --shift -100 --extsize 200 --broad;
  perl -lane 'print $F[0],"\t",$F[1],"\t",$F[2],"\t$ENV{thisn}\t",$F[4],"\t",$F[5]' /data/emarquez/peaks/${d}/$(basename ${f} .bed)_peaks.broadPeak  | bedtools sort > /data/emarquez/peaks/${d}/$(basename ${f} .bed)_peaks.bed;
done;

