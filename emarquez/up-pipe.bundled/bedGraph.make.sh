#!/bin/bash

# Submit to queue.
# This function creates a HOMER "Tag directory" in the directory tagdir in my data space. Tag directories are used to run various HOMER procedures; in this case they are used only to created bedGraph files.
# It creates a tag directory for each sample
# Output bedGraph files are placed in bedGraph directory
# Requires Homer

if [[ "$#" -lt 2 ]]
then
    echo "$(basename $0) [dirName] [groupName]"  1>&2
    echo "   [dirName]: Name of subdirectory where the BAM files are located" 1>&2
    echo "   [groupName]: Valid prefix of BAM file, depends on naming convention" 1>&2
    exit 1
fi

d=$(echo $1 | sed 's:/$::g');
g=$(echo $2 | sed 's:/$::g');
cd /data/emarquez/tagdir;

bamnames=$(find /data/emarquez/bamfiles/${d}/ -maxdepth 1 -name "${g}*.bam" | sort);
for f in ${bamnames}; do
	samp=$(basename $f | sed 's/\_/	/g' | perl -lane 'print $F[0],"_",$F[1]');
	sampname=$(basename $f .bam);
	[[ ! -d "/data/emarquez/tagdir/${d}/${samp}" ]] && mkdir "/data/emarquez/tagdir/${d}/${samp}";
	makeTagDirectory /data/emarquez/tagdir/${d}/${samp} ${f} -genome hg19 -single -fragLength 150;
	makeUCSCfile /data/emarquez/tagdir/${d}/${samp} -fragLength 150 > /data/emarquez/bedGraph/${d}/${sampname}.bedGraph;
done;
