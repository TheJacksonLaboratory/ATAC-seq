#!/bin/bash
#  Aligning fastq files to the genome and creating shifted_sorted.bam file -- preprocessing for ATAC-seq
#  ahead of peak calling
#  
#  STEP 4 of 6
#
#  Name   :  make_atac_seq_shifted_bam_3_shift_sam.sh
#
#  Purpose: This routine calls ATAC_BAM_shifter_gappedAlign.pl
#
#           ATAC_BAM_shifter_gappedAlign.pl  will take a BAM alignment file and adjust it as is appropriate for ATAC-seq:
#           
#           From Asli Uyar, PhD
#              "The bam file needs to be adjusted because Tn5 has a 9bp binding site, and it binds in the middle.  
#               Functionally that means that the DNA had to be accessible at least 4.5bp on either site of the insertion.
#
#               To my understanding the 4 for positive and 5 for negative was chosen at random, that could be wrong
#               for the negative strand the start position is the 5'-most, which doesn't actually change, it's the 
#               3'-position that changes, and to modify that you only need to change the read sequence and size.
#
#           If the read is on the positive strand (as determined by the sam flag) 
#
#                   1.  it will add 4 to the start, 
#                   2.  and subtract 5 from the partner start
#
#           if the read is on the negative strand, 
#
#                   1. 5 will subtracted from it's start 
#                   2. and 4 added to its mate start.
#
#            The length and read type will be adjusted in both cases and the read and quality string trimmed appropriately,
# 
# Example of how a BAM file will be altered
#
# Original:
# HWI-ST374:226:C24HPACXX:4:2214:15928:96004      99      chrI    427     255     101M    =       479     153 
# HWI-ST374:226:C24HPACXX:4:2214:15928:96004      147     chrI    479     255     101M    =       427     -153 
#
# Altered:
# HWI-ST374:226:C24HPACXX:4:2214:15928:96004      99      chrI    431     255     97M    =       474     144 
# HWI-ST374:226:C24HPACXX:4:2214:15928:96004      147     chrI    479     255     96M    =       431     -144 
#
#             
#  Author:
#            Asli Uyar, PhD
#            Shane Sanders, PhD
#            Anne Deslattes Mays, PhD
#
#  Date:     2016 October 10
#
#  Call:     make_atac_seq_shifted_bam_4_shift_sam.sh <ATAC-Seq Banchereau-Lab/GT-delivery/ATAC-seq directory (with trailing /> ]
#
#  Assumptions:  
#
#            1.  code has been checked out - fastqc has been run - working directory created
#            2.  we are in the same directory as the scripts.
#            3.  trimmomatic has been run
#            4.  bwa has been run
#
#-------------------------------------------------------------------------------------
#dataDIR=$1 #ARGV, contains folder with FASTQ files, right now needs trailing /

scriptDIR=$(pwd)
workingDIR=$scriptDIR/working
inputDIR=$workingDIR/trimmomatic/bwa
outputDIR=$workingDIR/trimmomatic/bwa

#find $workingDIR -name *R1_001.trim.fastq.gz > $workingDIR/filelist.txt
rm $workingDIR/samfilelist.txt
ls -1 $inputDIR/*.sam > $workingDIR/samfilelist.txt
FILENUMBER=$(wc -l $workingDIR/samfilelist.txt | cut -d' ' -f1)

#echo $FILENUMBER
rm $workingDIR/postBWA.qsub

echo \#!/bin/bash >> $workingDIR/postBWA.qsub
echo \#PBS -l nodes=1:ppn=16 >> $workingDIR/postBWA.qsub
echo \#PBS -l walltime=12:00:00 >> $workingDIR/postBWA.qsub
echo \#PBS -N shift_sam  >> $workingDIR/postBWA.qsub
echo \#PBS -t 1-$FILENUMBER >> $workingDIR/postBWA.qsub
echo module load python >> $workingDIR/postBWA.qsub
echo module load R >> $workingDIR/postBWA.qsub
echo module load perl/5.10.1 >> $workingDIR/postBWA.qsub
echo module load samtools/0.1.19 >> $workingDIR/postBWA.qsub
echo module load bedtools >> $workingDIR/postBWA.qsub
echo FILESAM=\$\(head -n \$PBS_ARRAYID $workingDIR/samfilelist.txt \| tail -1\) >> $workingDIR/postBWA.qsub
echo FILESORTED=\$\(basename "\${FILESAM}" \| sed \'s/\.sam/_sorted\.sam/g\'\) >> $workingDIR/postBWA.qsub
echo java -Xms1g -Xmx4g -jar /opt/compsci/picard/1.95/SortSam.jar INPUT=\$FILESAM OUTPUT=$outputDIR/\$FILESORTED SO=coordinate >>  $workingDIR/postBWA.qsub
echo FILERMDUP=\$\(basename \"\${FILESORTED}\" \| sed \'s/\.sam/_rmdup\.sam/g\'\) >> $workingDIR/postBWA.qsub
echo FILEMETRICS=\$\(basename \"\${FILESORTED}\" \| sed \'s/\.sam/_rmdup_metrics\.txt/g\'\) >>$workingDIR/postBWA.qsub
echo java -Xms1g -Xmx4g -jar /opt/compsci/picard/1.95/MarkDuplicates.jar INPUT=$outputDIR/\$FILESORTED OUTPUT=$outputDIR/\$FILERMDUP METRICS_FILE=$outputDIR/\$FILEMETRICS REMOVE_DUPLICATES=true >> $workingDIR/postBWA.qsub
echo FILEINSERT=\$\(basename \"\${FILERMDUP}\" \| sed \'s/\.sam/_insertSize\.txt/g\'\) >> $workingDIR/postBWA.qsub
echo FILEHISTO=\$\(basename \"\${FILERMDUP}\" \| sed \'s/\.sam/_insertSize\.pdf/g\'\) >> $workingDIR/postBWA.qsub
echo java -Xms1g -Xmx4g -jar /opt/compsci/picard/1.95/CollectInsertSizeMetrics.jar METRIC_ACCUMULATION_LEVEL=ALL_READS OUTPUT=$outputDIR/\$FILEINSERT HISTOGRAM_FILE=$outputDIR/\$FILEHISTO INPUT=$outputDIR/\$FILERMDUP >> $workingDIR/postBWA.qsub
echo FILENAME=\$\(basename \"\${FILERMDUP}\" \| sed \'s/\.sam/_shifted/g\'\) >> $workingDIR/postBWA.qsub
echo perl $scriptDIR/auyar/ATAC_BAM_shifter_gappedAlign.pl $outputDIR/\$FILERMDUP $outputDIR/\$FILENAME >> $workingDIR/postBWA.qsub

######
qsub -V $workingDIR/postBWA.qsub

echo "made the call to ATAC_BAM_shifter_gappedAlign.pl -- wait until that job is complete before going to STEP 3"
echo "use qstat -u <username> to check the status of your job" 
echo "done with STEP 4 when qsub is complete go to STEP 5 make_atac_seq_shifted_bam_5_bamtobed.sh"


