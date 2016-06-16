#!/bin/bash

# Submit this to queue.
# This function "scores" each BAM file in a given bamfiles directory on a given consensus peak file from the peaks directory.
# "Scoring" means counting raw reads for each BAM file on the regions defined by the consensus peak file.
# Requires bedtools
# Requires 3 input parameters, described below. These are outputs from rep.merge.sh, bam.index.sh, and peak.consensus.sh scripts

if [[ "$#" -lt 3 ]]
then
    echo "$(basename $0) [setName] [groupName] [minOver]"  1>&2
    echo "   [setName]: Name of subdirectory where the BAM files are located" 1>&2
    echo "   [groupName]: Valid prefix of BAM file, depends on naming convention" 1>&2
    echo "   [minOver]: minOverlap for consensus (e.g. 10)" 1>&2
    exit 1
fi

s=$(echo $1 | sed 's:/$::g');
g=$(echo $2 | sed 's:/$::g');
mover=$(echo $3 | sed 's:/$::g');
cd /data/emarquez/peaks/bamcounts;
bamnames=$(find /data/emarquez/bamfiles/${s}/ -maxdepth 1 -name "${g}*.bam" | sort);
> ${g}.bamnames1.tmp;
for f in ${bamnames};
 do
  basename $f | sed 's/\_/ /g' | perl -lane 'print $F[0],"_",$F[1]' >> ${g}.bamnames1.tmp;
done;
sed ':a;N;$!ba;s/\n/\t/g' ${g}.bamnames1.tmp | perl -ne 'print "chr\tstart\tend\tnhits\t",$_' > ${g}.bamnames2.tmp;
bedtools multicov -D -bams ${bamnames} -bed /data/emarquez/peaks/consensus/aging_sorted_merged_consensus_peaks_minOverlap${mover}.bed > \
		/data/emarquez/peaks/bamcounts/${g}_sorted_merged_peaks_rawcounts_minOVerlap${mover}.tmp;
cat ${g}.bamnames2.tmp /data/emarquez/peaks/bamcounts/${g}_sorted_merged_peaks_rawcounts_minOVerlap${mover}.tmp > \
        /data/emarquez/peaks/bamcounts/${g}_sorted_merged_peaks_rawcounts_minOVerlap${mover}.txt;
rm ${g}*.tmp;