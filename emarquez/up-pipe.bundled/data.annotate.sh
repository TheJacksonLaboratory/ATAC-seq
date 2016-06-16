#!/bin/bash

# This script reads a data matrix containing chr/start/end columns in position 1-3 and annotates these regions using chromHMM states and Homer annotations
# ChromHMM states downloaded from Roadmap epigenomics (this one is for PBMC)
# Requires Homer
# Data matrix also assumed to have a header
# Outputs three annotated files: _fullanno contains all Homer annotations, _annotated contains most commonly used fields; these two files have the same number of rows as the original data file
# The third output is a _chromHMM file. This file may have more rows than the data file, since multiple states can be assigned to a peak. In this case, peak coordinates are replicated as many times as states overlap them.

cd Dir1
chromHMM_dir=Dir2
data=data_matrix.txt

bed=$(basename $data .txt).bed # cuts columns 1-3 and removes header to build BED file
perl -ne 'print if $.>1' $data | cut -f1-3 > $bed
fullanno=$(basename $data .txt)_fullanno.txt
anno=$(basename $data .txt)_annotated

# chromHMM annotation
hmmanno_segments=$chromHMM_dir/E062_18_core_K27ac_segments.bed
head -1 $data | cut -f1-3 | sed 's/$/	chromHMMstate/g' > h.txt
bedtools intersect -wo -a $bed -b $hmmanno_segments | cut -f1-3,7 | sortBed | uniq > $(basename $data .txt)_chromHMM.tmp
sed s/E18/Quies/g $(basename $data .txt)_chromHMM.tmp \
	| sed s/E17/ReprPCWk/g | sed s/E16/ReprPC/g | sed s/E15/EnhBiv/g | sed s/E14/TssBiv/g | sed s/E13/Het/g | sed s/E12/ZNF_Rpts/g \
	| sed s/E11/EnhWk/g | sed s/E10/EnhA2/g | sed s/E9/EnhA1/g | sed s/E8/EnhG2/g | sed s/E7/EnhG1/g | sed s/E6/TxWk/g | sed s/E5/Tx/g \
	| sed s/E4/TssFlnkD/g | sed s/E3/TssFlnkU/g | sed s/E2/TssFlnk/g | sed s/E1/TssA/g | cat h.txt - \
	> $(basename $data .txt)_chromHMM.txt

# Homer annotation
annotatePeaks.pl $bed hg19 > $fullanno
cut -f2-4,8,10,16,19 $fullanno | sed 's/ //g' | sed s/\(.*\)//g | sed 's/\.[0-9]	/	/g' > ${anno}.tmp
head -1 ${anno}.tmp > h.txt
perl -ne 'print if $.>1' ${anno}.tmp | sortBed | perl -lane 'print $F[0],"\t",$F[1]-1,"\t",$F[2],"\t",$F[3],"\t",$F[4],"\t",$F[5],"\t",$F[6]' | cat h.txt - > ${anno}.txt
rm *tmp
rm h.txt