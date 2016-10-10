# make_atac_seq_shifted_bam
These scripts allow to allow the user to take ATAC-sequences through to shifted.bam.

There are 5 scripts that need to be run in succession -- each prior step must complete before the next.

Eventually we will separate out this part from the ATAC-seq repository -- but for now we shall leave it here
because of the dependencies that remain to the directory of scripts

Assumptions

# Workflow
      
     [ git clone the ATAC-seq repository ]
                 |
                 |
     [ make_atac_seq_shifted_bam_1_fastqc.sh <ATAC-Seq Banchereau-Lab/GT-delivery/ATAC-seq directory (with trailing /> ]
                 |
                 |
     [ make_atac_seq_shifted_bam_2_trimmomatic.sh <ATAC-Seq Banchereau-Lab/GT-delivery/ATAC-seq directory (with trailing /> ]
                 |
                 |
     [ make_atac_seq_shifted_bam_3_bwa.sh <ATAC-Seq Banchereau-Lab/GT-delivery/ATAC-seq directory (with trailing /> ]
                 |
                 |
     [ make_atac_seq_shifted_bam_4_shift_sam.sh <ATAC-Seq Banchereau-Lab/GT-delivery/ATAC-seq directory (with trailing /> ]
                 |
                 |
     [ make_atac_seq_shifted_bam_5_bamtobed.sh <ATAC-Seq Banchereau-Lab/GT-delivery/ATAC-seq directory (with trailing /> ]
                 |
                 |
     [ make_atac_seq_shifted_bam_6_cleanup.sh <ATAC-Seq Banchereau-Lab/GT-delivery/ATAC-seq directory (with trailing /> ]

#Get the scripts.
git clone this repository

  In the /data/Banchereau-Lab/GT-delivery/ATAC-seq directory, make a directory that contains a symbolic link to all the files
  that need to be processed 

#Call the first script

  make_atac_seq_shifted_bam_1_fastqc.sh
  Purpose:   
            
#Call the second script
  
  make_atac_seq_shifted_bam_2_trimmomatic.sh
  Purpose:  fastqc has been run now we run trimmomatic to identified the trimmend and untrimmed fastq files

  Call:     make_atac_seq_shifted_bam_1_run_fastqc.sh <ATAC-Seq Banchereau-Lab/GT-delivery/ATAC-seq directory (with trailing /> ]

  Assumptions:  

            1.  code has been checked out - fastqc has been run - working directory created
            2.  we are in the same directory as the scripts.


#Call the third script

  make_atac_seq_shifted_bam_3_bwa.sh 
  Purpose:  bwa is the alignment of the trimmed fastq files to the genome

  Call:     make_atac_seq_shifted_bam_1_run_fastqc.sh <ATAC-Seq Banchereau-Lab/GT-delivery/ATAC-seq directory (with trailing /> ]

  Assumptions:  

           1.  code has been checked out - fastqc has been run - working directory created
           2.  we are in the same directory as the scripts.
           3.  trimmomatic has been run

#Call the fourth script

  make_atac_seq_shifted_bam_4_shift_sam.sh

  Purpose: This routine calls ATAC_BAM_shifter_gappedAlign.pl

          
          From Asli Uyar, PhD
             "The bam file needs to be adjusted because Tn5 has a 9bp binding site, and it binds in the middle.  
              Functionally that means that the DNA had to be accessible at least 4.5bp on either site of the insertion.

              To my understanding the 4 for positive and 5 for negative was chosen at random, that could be wrong
              for the negative strand the start position is the 5'-most, which doesn't actually change, it's the 
              3'-position that changes, and to modify that you only need to change the read sequence and size.

          If the read is on the positive strand (as determined by the sam flag) 

                  1.  it will add 4 to the start, 
                  2.  and subtract 5 from the partner start

          if the read is on the negative strand, 

                  1. 5 will subtracted from it's start 
                  2. and 4 added to its mate start.

           The length and read type will be adjusted in both cases and the read and quality string trimmed appropriately,

 Example of how a BAM file will be altered

 Original:
 HWI-ST374:226:C24HPACXX:4:2214:15928:96004      99      chrI    427     255     101M    =       479     153 
 HWI-ST374:226:C24HPACXX:4:2214:15928:96004      147     chrI    479     255     101M    =       427     -153 

 Altered:
 HWI-ST374:226:C24HPACXX:4:2214:15928:96004      99      chrI    431     255     97M    =       474     144 
 HWI-ST374:226:C24HPACXX:4:2214:15928:96004      147     chrI    479     255     96M    =       431     -144 

 Call:     make_atac_seq_shifted_bam_4_shift_sam.sh <ATAC-Seq Banchereau-Lab/GT-delivery/ATAC-seq directory (with trailing /> ]

 Assumptions:  

            1.  code has been checked out - fastqc has been run - working directory created
            2.  we are in the same directory as the scripts.
            3.  trimmomatic has been run
            4.  bwa has been run


#Call the fifth script

  make_atac_seq_shifted_bam_5_bamtobed.sh

  Purpose: This routine sorts the shifted sam file and converts it to a bed file.

  Call:     make_atac_seq_shifted_bam_5_bamtobed.sh <ATAC-Seq Banchereau-Lab/GT-delivery/ATAC-seq directory (with trailing /> ]

  Assumptions:  

            1.  code has been checked out - fastqc has been run - working directory created
            2.  we are in the same directory as the scripts.
            3.  trimmomatic has been run
            4.  bwa has been run
            5.  shifted_bam has been called to make the sam file correct for ATAC_seq
           

#Call the sixth and final script -- this cleans up no longer needed files

  make_atac_seq_shifted_bam_6_cleanup.sh

  Purpose: This routine removes all the trim and trimU fastq files (intermediate files that may be reproduced)
           It also removes all sam files, and the unsorted bam file, also all .o and .e files for readibility.

 Call:     make_atac_seq_shifted_bam_6_cleanup.sh <ATAC-Seq Banchereau-Lab/GT-delivery/ATAC-seq directory (with trailing /> ]

 Assumptions:  

          1.  code has been checked out - fastqc has been run - working directory created
          2.  we are in the same directory as the scripts.
          3.  trimmomatic has been run
          4.  bwa has been run
          5.  shifted_bam has been called to make the sam file correct for ATAC_seq
          6.  the shifted_bam script and bedtobam has been called so all that is left to do is cleanup!

#What happens next?

  Now the real work begins -- peak calling, analysis, etc -- all the downstream work that is necessary
  to come to biologically meaningful conclusions.





