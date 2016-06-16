#!/bin/bash

while read line
do

dataPath1="$line"
cd $dataPath1
mkdir scripts
cp /home/auyar/ATAC_pipeline/mtDepletion/*.sh $dataPath1/scripts
cd scripts
qsub run_ATAC_pipeline.sh


done < sample.txt





