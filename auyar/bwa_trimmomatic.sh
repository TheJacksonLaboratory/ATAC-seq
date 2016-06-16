#!/bin/bash

dataPath=$1
ext1='trimmomatic/'
ext2='trimmomatic/bwa/'

inputPath=$dataPath$ext1
outputPath=$dataPath$ext2

for f in $(find "$inputPath" -name '*R1_001.trim.fastq')
do
	echo "Aligning: $f...\n"
	fileName2=$(basename "${f}" | sed 's/R1_001\.trim\.fastq/R2_001\.trim\.fastq/g')
	fileNameSAM=$(basename "${f}" | sed 's/R1_001\.trim\.fastq/hg19\.sam/g')
	fileNameMetrics=$(basename "${f}" | sed 's/R1_001\.trim\.fastq/metrics\.txt/g')
	/opt/compsci/bwa/0.7.12/bin/bwa mem -M /data/shared/genomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa $f $inputPath$fileName2 > $outputPath$fileNameSAM
done
echo "Finished mapping\n"

for f in $(find "$outputPath" -name '*.sam')
do
	echo "Shifting and cleaning $f then converting to bam... \n"
	fileSorted=$(basename "${f}" | sed 's/\.sam/_sorted\.sam/g')
	java -Xms1g -Xmx4g -jar /opt/compsci/picard/1.95/SortSam.jar INPUT=$f OUTPUT=$outputPath$fileSorted SO=coordinate
	fileRmdup=$(basename "${fileSorted}" | sed 's/\.sam/_rmdup\.sam/g')
	fileMetrics=$(basename "${fileSorted}" | sed 's/\.sam/_rmdup_metrics\.txt/g')
	java -Xms1g -Xmx4g -jar /opt/compsci/picard/1.95/MarkDuplicates.jar INPUT=$outputPath$fileSorted OUTPUT=$outputPath$fileRmdup METRICS_FILE=$outputPath$fileMetrics REMOVE_DUPLICATES=true
	# collect insert size info
	fileInsert=$(basename "${fileRmdup}" | sed 's/\.sam/_insertSize\.txt/g')
	fileHisto=$(basename "${fileRmdup}" | sed 's/\.sam/_insertSize\.pdf/g')
	java -Xms1g -Xmx4g -jar /opt/compsci/picard/1.95/CollectInsertSizeMetrics.jar METRIC_ACCUMULATION_LEVEL=ALL_READS OUTPUT=$outputPath$fileInsert HISTOGRAM_FILE=$outputPath$fileHisto INPUT=$outputPath$fileRmdup
	# now shift and filter for mapping quality
	fileName=$(basename "${fileRmdup}" | sed 's/\.sam/_shifted/g')
	perl /home/auyar/scripts/ATAC_BAM_shifter_gappedAlign.pl $outputPath$fileRmdup $outputPath$fileName

done
echo "Finished shifting, cleaning and converting\n"


for f in $(find "$outputPath" -name '*shifted.bam')
do
	echo "process\n";
	fileName=$(basename "${f}" | sed 's/\.bam/_sorted.bam/g')
	fileName1=$(basename "${f}" | sed 's/\.bam/_sorted/g')
	fileName2=$(basename "${f}" | sed 's/\.bam/_sorted.bed/g')
	echo "\n Sorting and Indexing file: $f... \n"
	samtools sort $f $outputPath$fileName1
	samtools index $outputPath$fileName
	bedtools bamtobed -i $outputPath$fileName > $outputPath$fileName2
done

