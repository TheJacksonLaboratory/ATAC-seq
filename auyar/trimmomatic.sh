#!/bin/bash

dataPath=$1
ext='trimmomatic/'
inputPath=$dataPath$ext
outputPath=$dataPath$ext

for f in $(find "$dataPath" -name '*R1_001.fastq')
do
	echo ": $f...\n"
	fileName2=$(basename "${f}" | sed 's/R1_001\.fastq/R2_001\.fastq/g')
	echo $f
	echo $fileName2
	fileName3=$(basename "${f}" | sed 's/R1_001\.fastq/R1_001\.fastq_filtered/g')
	fileName4=$(basename "${fileName2}" | sed 's/R2_001\.fastq/R2_001\.fastq_filtered/g')
	fileName5=$(basename "${f}" | sed 's/R1_001\.fastq/R1_001\.trimU\.fastq/g')
	fileName6=$(basename "${fileName2}" | sed 's/R2_001\.fastq/R2_001\.trimU\.fastq/g')
	java -jar /opt/compsci/Trimmomatic/0.33/trimmomatic-0.33.jar PE $f $dataPath$fileName2 $outputPath$fileName3 $outputPath$fileName5 $outputPath$fileName4 $outputPath$fileName6 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 
done
echo "Finished QT & QF\n"


for f in $(find "$inputPath" -name '*R1_001.fastq_filtered')
do
	echo "Adapter trimming: $f...\n"
	fileName2=$(basename "${f}" | sed 's/R1_001\.fastq_filtered/R2_001\.fastq_filtered/g')
	python /home/auyar/scripts/pyadapter_trim.py -a $f -b $inputPath$fileName2
done
echo "Finished adapter trimming\n"

