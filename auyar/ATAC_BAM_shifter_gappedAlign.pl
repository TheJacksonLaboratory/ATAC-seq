#!/usr/bin/perl
use strict;
use warnings;


# This program will take a BAM alignment file and adjust it as is appropriate for ATAC-seq:
# If the read is on the positive strand (as determined by the sam flag) it will add 4 to the start, and subtract 5 from the partner start
# if the read is on the negative strand, 5 will subtracted from it's start and 4 added to its mate start.
# the length and read type will be adjusted in both cases and the read and quality string trimmed appropriately

# this is done because Tn5 has a 9bp binding site, and it binds in the middle.  Functionally that means that the DNA had to be accessible at least 4.5bp on either site of the insertion.
# to my understanding the 4 for positive and 5 for negative was chosen at random, that could be wrong
# for the negative strand the start position is the 5'-most, which doesn't actually change, it's the 3'-position that changes, and to modify that you only need to change the read sequence and size.

# example of how a BAM file will be altered
# Original:
# HWI-ST374:226:C24HPACXX:4:2214:15928:96004      99      chrI    427     255     101M    =       479     153 # should also have read and quality, but were omitted for readability
# HWI-ST374:226:C24HPACXX:4:2214:15928:96004      147     chrI    479     255     101M    =       427     -153 # should also have read and quality, but were omitted for readability
#
# Altered:
# HWI-ST374:226:C24HPACXX:4:2214:15928:96004      99      chrI    431     255     97M    =       474     144 # read and quality were moitted for readability, but they would each be trimmed 4bp
# HWI-ST374:226:C24HPACXX:4:2214:15928:96004      147     chrI    479     255     96M    =       431     -144 # read and quality were moitted for readability, but they would each be trimmed 5bp


# make sure the right arguments are being passed
die "USEAGE:
This program will accept alignment files from an ATAC-seq experiment and adjust the mapping fragments to better reflect the biology, outputting a modified bam file
Arguments required: 1. sam or bam alignment file (must end in .bam or .sam) 2. Filehandle for the output
" unless(($ARGV[0] =~ m/\.sam$/ ||$ARGV[0] =~ m/\.bam$/) && $ARGV[1] =~ m/\w/);


my $bamFile = shift;
my $out = shift;
my $outBam = $out.".bam";
my $temp = $out.".sam";

# only reads mapping with a quality of 10 will be considered
if ($bamFile=~m/.bam$/){
    die "cannot open $bamFile.\n" unless(open BAM,"samtools view -h -q 10 $bamFile |");
}elsif($bamFile=~m/.sam$/){
    die "cannot open $bamFile.\n" unless(open BAM,"samtools view -h -S -q 10 $bamFile |");
}else{
    die "Input alignment file is not a .sam or .bam file\n";
}

die "cannot open $temp.\n" unless(open TEMP,">$temp");

my %flags=samFlags(); # this is a little sub that generates a hash with keys of 'unmapped', 'positive' and 'negative'; that is, if the read maps, and if so, to which strand, regardless of anything else.

while(<BAM>){
    chomp;  # removing new line
    # adding the header lines (if you used -h in the samools command), we'll want it later, and going on to the reads
    print TEMP $_."\n" if(/^(\@)/);
    next if(/^(\@)/);

    my @line = split(/\t/);  # splitting SAM line into array

    # skip reads if both both pairs were not
    # the specific read is unmapped,
    # or map to the mitochondria.
    next unless ($line[1] ~~ $flags{'properlyPaired'});
    #next if($line[2] ~~ m/^chrM/ || $line[2] ~~ m/_random$/ || $line[2] ~~m/^chrUn/);
    next if($line[2] ~~ m/_random$/ || $line[2] ~~m/^chrUn/);
    if($line[8]==0){ # just in case
	print STDERR $line[0]." has a length of 0\n$_\n"; # just in case
    }

    if($line[1] ~~ $flags{'negative'}){ # if the read is mapped to the reverse strand
	$line[7]=$line[7]+4; # move the mate start 4 bp
	$line[8]=$line[8]+9; # adjust the inferred insert size
	$line[9] =substr($line[9],0,-5); # remove 5 bp from the reads
	$line[10] =substr($line[10],0,-5);# remove the 5bp from quality as well

	# this is ugly, but this is a way to change the CIGAR to reflect that I'm shortening the read
	# for the negative stran I want to subtract the 5p from the last entry (which should be M, for match)
	$line[5] = &CIGARtrim($line[5],'-',5); # from the CIGAR string, trim 5, and we're on the negtive strand

    }elsif($line[1] ~~ $flags{'positive'}){ # if the read is mapped to the positive strand
	$line[3]=$line[3]+4; # move the start 4 bp
	# these were removed, because they were incorrect, you don't actually move this end of the negative read (which this would be, because it's the mate of a positive strand read).
	# for a negative read you just trim from the other end, the start position stays the same.
	#$line[7]=$line[7]-5; # move the mate start 5 bp 
	#$line[7]=1 if $line[7]<=0; # just to make sure we didn't go off the end
	$line[8]=$line[8]-9; # adjust the inferred insert size
	$line[9] =substr($line[9],4);# remove 4bp from the reads
	$line[10] =substr($line[10],4);# remove the 4bp from quality as well

	# this is ugly, but this is a way to change the CIGAR to reflect that I'm shortening the read
	# for the positive strand I want to subtract the 4 from the first entry (which should be M, for match)
	$line[5] = &CIGARtrim($line[5],'+',4); # from the CIGAR string, trim 4, and we're on the positive strand

    }

    print TEMP join("\t",@line)."\n";
}

close BAM;
close TEMP;

# convert the sam file to the bam file I want, and clean up
`samtools view -Sb $temp > $outBam`;
`rm $temp`;

exit;


##### Sub routines #####

# this generates a hash filled with the most common sam flags and wether the reads are mapped, and if so to which strand
# NOTE this does not address any of the other issues that can be interpretted from the flag (correct orientation within or ourside of the insert size)
sub samFlags{
    my %flags=(
	'unmapped' => [133,165,181,101,117,69,77,141],
	'properlyPaired'=>[99,147,83,163,81,161,97,145],
	'improperInsert' => [81,161,97,145],
	'positive' => [73,89,137,99,163,67,131,161,97,65,129],
	'negative' => [121,153,185,147,83,115,179,81,145,113,177]
	);
    return (%flags);
}

# this cuts off the first number and letter from a CIGAR string
# and returns the trimmed CIGAR plus the number and letter trimmed off
sub trimFirstCIGARunit{
    my $cigar = shift;
    my $firstChar = (split/\d+/, $cigar)[1];
    #print STDERR "$cigar\n$firstChar\n";
    my $firstNum = (split/$firstChar/, $cigar)[0];
    my $length =length($firstNum.$firstChar);
    return(substr($cigar,$length), $firstNum, $firstChar);
}

# this cuts off the last number and letter from a CIGAR string
# and returns the trimmed CIGAR plus the number and letter trimmed off
sub trimLastCIGARunit{
    my $cigar = shift;
    my $lastChar = substr($cigar,-1); # split was acting weird, so this should grab the last character
    my $lastNum = (split/[MINDS]/, $cigar)[-1]; # these are the possible letters equalling Match, Insert, N, Deletion, Shoft clip
    my $length =length($lastNum.$lastChar);
    return(substr($cigar,0,(-1*$length)), $lastNum, $lastChar);
}

# this will take a cigar string, the strand the read maps to and a number of reads to cut off
# it returns the properly trimmed CIGAR string
sub CIGARtrim{
    my ($cigar, $strand, $num) = @_;
    my $toReturn;
	if ($strand ~~ '-'){
	    my ($cigarToKeep, $lastNum, $lastChar) = &trimLastCIGARunit($cigar);
	    my $newLastNum = $lastNum - $num;

	    if ($lastChar ~~ 'D'){ # skip to the next letter because the deletion is completely skipped over
		($cigarToKeep, $lastNum, $lastChar) = &trimLastCIGARunit($cigarToKeep);
		$newLastNum = $lastNum - $num;
	    }

	    if ($newLastNum==0){
		my $nextChar = substr($cigarToKeep,-1);
		if ($nextChar ~~ 'D'){ # skip to the next letter because a read shouldn't end in a deletion
		    ($cigarToKeep, $lastNum, $lastChar) = &trimLastCIGARunit($cigarToKeep);
		}
		$toReturn = $cigarToKeep;
	    }elsif($newLastNum<0){ # Not everything has been trimmed off, so iterate over this again
		$toReturn=&CIGARtrim($cigarToKeep,$strand,-1*$newLastNum);
	    }else{
		$toReturn=$cigarToKeep.$newLastNum.$lastChar;
	    }

	}elsif($strand ~~ '+'){
	   # print STDERR "$cigar\n";
	    my ($cigarToKeep, $firstNum, $firstChar) = &trimFirstCIGARunit($cigar);
	    my $newFirstNum = $firstNum - $num;

	    if ($firstChar ~~ 'D'){ # skip to the next letter because the deletion is completely skipped over
		($cigarToKeep, $firstNum, $firstChar) = &trimFirstCIGARunit($cigarToKeep);
		$newFirstNum = $firstNum - $num;
	    }

	    if ($newFirstNum==0){
		my $nextChar = substr($cigarToKeep,-1);
		if ($nextChar ~~ 'D'){ # skip to the next letter because a read shouldn't end in a deletion
		    ($cigarToKeep, $firstNum, $firstChar) = &trimFirstCIGARunit($cigarToKeep);
		}
		$toReturn = $cigarToKeep;
	    }elsif($newFirstNum<0){ # Not everything has been trimmed off, so iterate over this again
		$toReturn=&CIGARtrim($cigarToKeep,$strand,-1*$newFirstNum);
	    }else{
		$toReturn=$newFirstNum.$firstChar.$cigarToKeep;
	    }

	}else{
	    die ("strand in CIGARtrim must be '+' or '-'\n");
	}
	return ($toReturn);
}
