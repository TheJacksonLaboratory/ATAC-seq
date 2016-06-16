#!/bin/bash

# Submit from queue
# This function uses PIQ to estimate footprint locations in BAM files
# syntax: footprints_PBMC_full.sh bamFile oDir fpStart fpEnd

if [[ "$#" -lt 4 ]]
then
    echo "$(basename $0) [bamFile] [oDir]"  1>&2
    echo "   [bamFile]: path to raw read file (BAM format)" 1>&2
    echo "   [oDir]: output directory" 1>&2
    echo "   [fpStart]: first FP to read from the FP list" 1>&2
    echo "   [fpEnd]: last FP to read from the FP list" 1>&2
    exit 1
fi

bamFile=$(echo $1 | sed 's:/$::g')
oDir=$(echo $2 | sed 's:/$::g')
fp0=$(echo $3)
fp1=$(echo $4)

fileName=$(basename "${bamFile}" |  sed 's/\.bam//g');
fileNamePIQPOOLED=${oDir}/outputFiles/$(basename ${fileName} | sed 's/\.bed//g')_pooledPIQcalls.bed
dirNamePIQ="${oDir}/$(basename ${oDir}_piq)";
dirNamePIQCALLS="${oDir}/$(basename ${oDir}_piqcalls)";

[[ ! -d "${oDir}" ]] && mkdir "${oDir}"
[[ ! -d "${oDir}/outputFiles" ]] && mkdir "${oDir}/outputFiles"
[[ ! -d "${oDir}/outputFiles/PIQ-single" ]] && mkdir "${oDir}/outputFiles/PIQ-single"
[[ ! -d "${dirNamePIQ}" ]] && mkdir "${dirNamePIQ}";
[[ ! -d "${dirNamePIQCALLS}" ]] && mkdir "${dirNamePIQCALLS}";

## PIQ-specific directories
## Location of the common.r file that does package load and defines parameters.
commonfile="~/seqanal/PIQ/common.r";
## Directory in which the PWM files were outputted using pwmmatch
pwmdir="~/seqanal/PIQ/motif.matches/hg19.masked/"
## Temporary read/write dir. Should be fast and local to avoid bottlenecking
tmpdir="${oDir}/tmp/"
[[ ! -d "${tmpdir}" ]] && mkdir "${tmpdir}";
## Processed bamfile created by bam2rdata script.
rbamfile="${dirNamePIQ}/${fileName}.RData"

## Loop PIQ over TF motifs to estimate putative footprints
Rscript ~/seqanal/PIQ/bam2rdata.r $commonfile $rbamfile $bamFile
for pwmid in `seq ${fp0} ${fp1}`
do
	Rscript ~/seqanal/PIQ/pertf.r $commonfile $pwmdir $tmpdir $dirNamePIQ $rbamfile $pwmid
	##clean up
	rm ${tmpdir}/${pwmid}/*
	rmdir ${tmpdir}/${pwmid}
	rm ${tmpdir}/*
done

## Clean up previous runs
rm ${oDir}/outputFiles/*
rm ${oDir}/outputFiles/PIQ-single/*
rm ${fileNamePIQPOOLED}.tmp
rm ${fileNamePIQPOOLED}

##Sorts by chr>start output from PIQ (all calls), filtering for purity >= 0.9, saving as *calls.bed file

for f in $(find "${dirNamePIQ}" -name '*-calls.all.bed'  -maxdepth 1 | grep -v '.RC')
do
	topcallsFile=${dirNamePIQCALLS}/$(basename "${f}" |  sed 's/\-calls\.all\.bed/-calls\.bed/g');
	allcallsRCFile=${dirNamePIQ}/$(basename "${f}" |  sed 's/\-calls\.all\.bed/\.RC\-calls\.all\.bed/g');
	export ucallsname=$(basename ${topcallsFile} | sed 's/\.bed//g' | sed 's/\-/\_/g');
 	export ucallstag=$(basename ${topcallsFile} | sed 's/\-calls.bed//g' | sed 's/\-/\_/g');
	cat ${f} ${allcallsRCFile} | bedtools sort | perl -lane 'print $_ if($F[4]>=900) ' | perl -ple 'print "track name=$ENV{ucallsname} description=\" Top PIQ calls for $ENV{ucallstag} \" useScore=1" if $. == 1' > ${topcallsFile}
	cat ${topcallsFile} >> ${fileNamePIQPOOLED}.tmp;
#### Note: the above step can remove all entries from the calls.bed file, in which case the following lines will return an error and the corresponding outputs will be empty;
## 	echo "$f"
done

## removes track lines from file
perl -ne 'print unless /^t/' ${fileNamePIQPOOLED}.tmp > ${fileNamePIQPOOLED}
rm ${fileNamePIQPOOLED}.tmp

## Copy non-empty individual TF PIQ calls to output singles folder
for g in $(find "${dirNamePIQCALLS}" -name '*bed'  -maxdepth 1);
do
	fsize=$(du $g | awk '{print $1}');
 	[[ $fsize != 0 ]] && cp $g ${oDir}/outputFiles/PIQ-single;
done

## Compress to go (kept in base sample directory)
tar -cf ${oDir}/outputFiles_$(basename "${fileName}" |  sed 's/$/\.tar/g') ${oDir}/outputFiles/*
gzip ${oDir}/outputFiles_$(basename "${fileName}" |  sed 's/$/\.tar/g')