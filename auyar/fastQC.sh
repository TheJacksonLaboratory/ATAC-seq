#!/bin/bash

dataPath=$1
ext='fastQC/'
outputPath=$dataPath$ext


for f in $(find "$dataPath" -name '*R1_001.fastq')
do
	echo "QC: $f...\n"
	fileName2=$(basename "${f}" | sed 's/R1_001\.fastq/R2_001\.fastq/g')
	echo $f
	echo $fileName2
	/opt/compsci/FastQC/0.11.3/fastqc $f $dataPath$fileName2 -o $outputPath 
done
echo "Finished QC\n"


