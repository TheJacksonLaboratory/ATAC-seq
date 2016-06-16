#!/bin/bash

dataPath=$1
ext='trimmomatic/bwa/'
ext2='trimmomatic/bwa/macs2/'
outputPath=$dataPath$ext2


for f in $(find "$dataPath$ext" -name '*_sorted.bed')
do
	echo "process\n";
	echo "f: $f \n"
	fileName=$(basename "${f}" | sed 's/\.bed//g')
	echo "fileName: $fileName \n "
	/opt/compsci/MACS/2.1.0.20151222/bin/macs2 callpeak -t $f -f BED -n $inputPath$fileName -g 'hs' --nomodel --shift -100 --extsize 200 -B --broad --outdir $outputPath
done


