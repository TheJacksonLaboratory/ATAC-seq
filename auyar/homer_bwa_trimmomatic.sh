#!/bin/bash


dataPath=$1
inputExt='trimmomatic/bwa/'
inputPath=$dataPath$inputExt
outputExt='trimmomatic/bwa/homer/'
outputPath=$dataPath$outputExt

for f in $(find "$inputPath" -name '*sorted.bed')
do
	fileName3=$(basename "${f}" | sed 's/\.bed/\_TAG/g')
	#fileName3=$fileName2"_TAG"
	fileName4=$(basename "${f}" | sed 's/\.bed/\.bedGraph/g')
	echo "$f, $fileName3, $fileName4"
	/opt/compsci/homer/4.6/bin/makeTagDirectory $outputPath$fileName3  $f -format bed -single -genome hg19
	/opt/compsci/homer/4.6/bin/makeUCSCfile $outputPath$fileName3 > $outputPath$fileName4
done
echo "Finished homer...\n"
