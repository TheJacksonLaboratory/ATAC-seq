#!/bin/bash

# Submit from queue
# This function uses Homer and Homer tag directory data to build footprint histograms for footprint locations found by PIQ
# This function runs for a single sample (multiple footprints, single tag directory)
# Clips BAM reads to 1bp
# syntax: footprints_PBMC_full.sh bamFile oDir hDir
# Current hDir in my data space is 'histograms'
# Creates a histogram directory for each sample
# Requires Homer annotatePeaks.pl to build histogram
# Depends on output from bedGraph.make.sh (tag directory) and footprint.call.sh (PIQ footprints)

if [[ "$#" -lt 3 ]]
then
    echo "$(basename $0) [bamFile] [oDir]"  1>&2
    echo "   [bamFile]: path to raw read file (BAM)" 1>&2
    echo "   [oDir]: output directory" 1>&2
    echo "   [hDir]: histogram directory" 1>&2
    exit 1
fi

bamFile=$(echo $1 | sed 's:/$::g')
oDir=$(echo $2 | sed 's:/$::g')
hDir=$(echo $3 | sed 's:/$::g')

fileName=$(basename "${bamFile}" |  sed 's/\.bam//g');
dirNamePIQCALLS="${oDir}/$(basename ${oDir}_piqcalls)";
dirNamePIQHIST="${hDir}/Dir1"; # Dir1=group of samples
dirNameTAG="/data/emarquez/tagdir/Dir1"; # Dir1=group of samples
let flank=1000; # flank define the size in bp of the region at each side of a footprint that will be sampled to build the histogram

##builds a 1bp resolution read histogram centered at footprints

## USE FOR SPECIFIC TF MOTIFS:
# let flank=1000;
# for i in FOX PBX PAX CTCF NFYA HNF MAF NKX
# do
# 	for f in $(find "${dirNamePIQCALLS}" -name "*${i}*-calls.bed"  -maxdepth 1)
# 	do
#     	echo "$f"
# 		topcallsHist=${dirNamePIQHIST}/$(basename "${f}" |  sed 's/\-calls\.bed/\-hist\.txt/g');
# 		echo "Building histogram for $f"
#     	annotatePeaks.pl ${f} hg19 -size $flank -hist 1 -ghist -d $dirNameTAG > $topcallsHist;
# 	done
# done

## USE FOR WHOLESALE TF MOTIFS:
for f in $(find "${dirNamePIQCALLS}" -name "*-calls.bed"  -maxdepth 1)
do
   	echo "$f"
	topcallsHist=${dirNamePIQHIST}/$(basename "${f}" |  sed 's/\-calls\.bed/\-hist\.txt/g');
	echo "Building histogram for $f"
   	annotatePeaks.pl ${f} hg19 -size $flank -hist 1 -ghist -d $dirNameTAG > $topcallsHist;
done