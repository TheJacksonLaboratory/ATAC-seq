dataPath=$1
ext1='trimmomatic/'
ext2='trimmomatic/bwa/'

scriptDIR=$(pwd)
inputPath=$dataPath$ext1
outputPath=$dataPath$ext2

#!/bin/bash

dataDIR=$1 #ARGV, contains folder with FASTQ files, right now needs full path
#find $dataDIR -name *R1_001.trim.fastq.gz > $dataDIR/filelist.txt
rm $dataDir/samfilelist.txt
ls -1 $dataDIR$ext2 > $dataDIR/samfilelist.txt
FILENUMBER=$(wc -l $dataDIR/samfilelist.txt | cut -d' ' -f1)

#echo $FILENUMBER
rm $dataDIR/postBWA.qsub

echo \#!/bin/bash >> $dataDIR/postBWA.qsub
echo \#PBS -l nodes=1:ppn=16 >> $dataDIR/postBWA.qsub
echo \#PBS -l walltime=48:00:00 >> $dataDIR/postBWA.qsub
echo \#PBS -N bwa  >> $dataDIR/postBWA.qsub
echo \#PBS -t 1-$FILENUMBER >> $dataDIR/postBWA.qsub
echo module load python >> $dataDIR/postBWA.qsub
echo module load R >> $dataDIR/postBWA.qsub
echo module load perl/5.10.1 >> $dataDIR/postBWA.qsub
echo module load samtools/0.1.19 >> $dataDIR/postBWA.qsub
echo module load bedtools >> $dataDIR/postBWA.qsub
echo FILESAM=\$\(head -n \$PBS_ARRAYID $outputDIR/samfilelist.txt \| tail -1\) >> $dataDIR/postBWA.qsub
echo FILESORTED=\$\(basename "\${FILESAM}" \| sed \'s/\.sam/_sorted\.sam/g\'\) >> $dataDIR/postBWA.qsub
echo java -Xms1g -Xmx4g -jar /opt/compsci/picard/1.95/SortSam.jar INPUT=\$FILESAM OUTPUT=$outputPath\$FILESORTED SO=coordinate >>  $dataDIR/postBWA.qsub
echo FILERMDUP=\$\(basename \"\${FILESORTED}\" \| sed \'s/\.sam/_rmdup\.sam/g\'\) >> $dataDIR/postBWA.qsub
echo FILEMETRICS=\$\(basename \"\${FILESORTED}\" \| sed \'s/\.sam/_rmdup_metrics\.txt/g\'\) >>$dataDIR/postBWA.qsub
echo java -Xms1g -Xmx4g -jar /opt/compsci/picard/1.95/MarkDuplicates.jar INPUT=$outputPath\$FILESORTED OUTPUT=$outputPath\$FILERMDUP METRICS_FILE=$outputPath\$FILEMETRICS REMOVE_DUPLICATES=true >> $dataDIR/postBWA.qsub
echo FILEINSERT=\$\(basename \"\${FILERMDUP}\" \| sed \'s/\.sam/_insertSize\.txt/g\'\) >> $dataDIR/postBWA.qsub
echo FILEHISTO=\$\(basename \"\${FILERMDUP}\" \| sed \'s/\.sam/_insertSize\.pdf/g\'\) >> $dataDIR/postBWA.qsub
echo java -Xms1g -Xmx4g -jar /opt/compsci/picard/1.95/CollectInsertSizeMetrics.jar METRIC_ACCUMULATION_LEVEL=ALL_READS OUTPUT=$outputPath$FILEINSERT HISTOGRAM_FILE=$outputPath\$FILEHISTO INPUT=$outputPath\$FILERMDUP >> $dataDIR/postBWA.qsub
echo FILENAME=\$\(basename \"\${FILERMDUP}\" \| sed \'s/\.sam/_shifted/g\'\) >> $dataDIR/postBWA.qsub
echo perl $scriptDIR/auyar/ATAC_BAM_shifter_gappedAlign.pl $outputPath\$FILERMDUP $outputPath\$FILENAME >> $dataDIR/postBWA.qsub

##
##
#
#for f in $(find "$outputPath" -name '*shifted.bam')
#do#
#	echo "process\n";
#	fileName=$(basename "${f}" | sed 's/\.bam/_sorted.bam/g')
#	fileName1=$(basename "${f}" | sed 's/\.bam/_sorted/g')
#	fileName2=$(basename "${f}" | sed 's/\.bam/_sorted.bed/g')
#	echo "\n Sorting and Indexing file: $f... \n"
#	samtools sort $f $outputPath$fileName1
#	samtools index $outputPath$fileName
#	bedtools bamtobed -i $outputPath$fileName > $outputPath$fileName2
#done

