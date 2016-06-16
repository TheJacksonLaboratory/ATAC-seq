#!/bin/bash

cd ..
gunzip *R1_001.fastq.gz
nRawReads4=`wc -l < *R1_001.fastq`
nRawReads=$((nRawReads4/4))

cd trimmomatic
gunzip *R1_001.fastq_filtered.gz
nQualityFilteredReads4=`wc -l < *R1_001.fastq_filtered`
nQualityFilteredReads=$((nQualityFilteredReads4/4))
gunzip *R1_001.trim.fastq.gz
nAdapterFilteredReads4=`wc -l < *R1_001.trim.fastq`
nAdapterFilteredReads=$((nAdapterFilteredReads4/4))

#cat *R1_ALL.fastq_filtered | awk '{if(NR%4==2) print length($1)}' > $dataPath/trimmomatic/readLengthQT.txt
#cat *R1_ALL.trim.fastq | awk '{if(NR%4==2) print length($1)}' > $dataPath/trimmomatic/readLengthAT.txt

#/opt/compsci/R/3.1.1/bin/Rscript /home/auyar/ATAC_pipeline/dataStats/plotReadLengthDistn.r $dataPath/trimmomatic/

#gzip *R1_ALL.fastq_filtered
#gzip *R1_ALL.trim.fastq

cd bwa
percentDuplication=`head -8 *metrics.txt | tail -1 | cut -f8`
nAlignedReads2=`wc -l < *.bed`
nAlignedReads=$((nAlignedReads2/2))

alignmentRate=$(echo "scale=4; $nAlignedReads/$nRawReads" | bc)

cd macs2
nPeaks=`wc -l < *broadPeak`

cd ../../../

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "RawReads" "QualityFilteredReads" "AdapterFilteredReads" "PercentDuplication" "AlignedReads" "AlignmentRate" "Peaks" >> nReads.txt  

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$nRawReads" "$nQualityFilteredReads" "$nAdapterFilteredReads" "$percentDuplication" "$nAlignedReads" "$alignmentRate" "$nPeaks" >> nReads.txt



