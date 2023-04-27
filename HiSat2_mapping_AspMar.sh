#!/bin/bash

#This script is made by bash and used for RNA-seq. 
#It first go through the fastqc. 

#load modules
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load sra
module load fastqc/0.10.1
module load trimmomatic/0.39
module load hisat2
module load stringtie/2.2.1
module load gffcompare
module load python/2.7.1
module load gcc/9.3.0
module load samtools
module load bcftools/1.2
module load gffread/
module load gffcompare/
module load multiqc/
module load trim_galore/

#record the start date for calculating how long we used. 
start=`date +%s`
echo $start

#please provide your user name 
MyID=$(cat ./input.txt)
echo $MyID 

#create the txt file in the directory which include the name of the reference species
#REF=$(cat ./reference_species.txt)

REF_origin=AspMar
REF_other=lizard

#make directory
DD=/scratch/$MyID/FunGenShare2023/WhiptailRawData
WD=/scratch/$MyID/PROJECT2
CD=/scratch/$MyID/PROJECT/CleanData
#CD2=/scratch/$MyID/PROJECT/trimgalore_Cleandata_after
REFD_origin=/scratch/$MyID/PROJECT2/$REF_origin
REFD_other=/scratch/$MyID/PROJECT2/$REF_other
MAPD_other=/scratch/$MyID/PROJECT2/Map_other
MAPD_other_nonannotation=/scratch/$MyID/PROJECT2/Map_other_nonannotation
MAPD_origin=/scratch/$MyID/PROJECT2/Map_origin
COUNTSD_other=/scratch/$MyID/PROJECT2/Counts_StringTie_other
COUNTSD_other_nonannotation=/scratch/$MyID/PROJECT2/Counts_StringTie_other_nonannotation
COUNTSD_origin=/scratch/$MyID/PROJECT2/Counts_StringTie_origin
RESULT_other=/scratch/$MyID/PROJECT2/Counts_other
RESULT_other_nonannotation=/scratch/$MyID/PROJECT2/Counts_other_nonannotation
RESULT_origin=/scratch/$MyID/PROJECT2/Counts_origin
#CS=/scratch/$MyID/PROJECT/PreCleanQuality
#CS2=/scratch/$MyID/PROJECT/PostCleanQuality
#CS22=/scratch/$MyID/PROJECT/trimgalore_PostCleanQuality

#-p tells it to make any upper level directories that are not there, notice how this will also make the WD 
mkdir $WD
#mkdir $CD
mkdir $REFD_origin
mkdir $REFD_other

mkdir $MAPD_origin
mkdir $MAPD_other
mkdir $MAPD_other_nonannotation

mkdir $COUNTSD_origin
mkdir $COUNTSD_other
mkdir $COUNTSD_other_nonannotation

mkdir $RESULT_origin
mkdir $RESULT_other
mkdir $RESULT_other_nonannotation

#mkdir $CS
#mkdir $CS2
#mkdir $CD2
#mkdir $CS22
#move to the Data Directory 
cd $DD 

#######download data files from NCBI: SRA using the RUN IDs 
##from SRA use the SRA tool kit See NCBI website https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/
##this downloads the SRA file and converts to fastq
##-F Defline contains only original sequence name. 
##-I Append read id after spot id as 'accession.spot.readid' on defline.
##splits the files into R1 and R2 (forward reads, reverse reads)

#vdb-config --interactive 
#THIS COMMAND IS USED FOR CONNECT WITH SUPTERCOMPUTER, I AM NOT SURE WE NEED IT OR NOT, MAYBE ONLY USED FOR THE FIRST TIME WE CONNECT 
#fastq-dump -F --split-files SRR6819016
#fastq-dump -F --split-files SRR6819022


##### Extra ####
## If you are downloaded data from a sequencing company instead of NCBI, using wget for example, then calculate the md5sum values of all the files in the folder (./*), and read into a text file.
## then you can compare the values in this file with the ones provided by the company.
#md5sum ./* > md5sum.txt

##### Extra ####
## If you data comes with multiple R1 and R2 files per individual. You can contatenate them together using "cat" before running FASTQC
## see examples below for one file. You will probably want to use a loop to process through all the files.
#cat SRR6819014*_R1_*.fastq.gz > SRR6819014_All_R1.fastq.gz
#cat SRR6819014*_R2_*.fastq.gz > SRR6819014_All_R2.fastq.gz


############## FASTQC to assess quality of the sequence data
## FastQC: run on each of the data files that have 'All' to check the quality of the data
## The output from this analysis is a folder of results and a zipped file of results


############## FASTQC to assess quality of the sequence data
## FastQC: run on each of the data files that have 'All' to check the quality of the data
## The output from this analysis is a folder of results and a zipped file of results
#fastqc * --outdir=$CS


#######  Tarball the directory containing the FASTQC results so we can easily bring it back to our computer to evaluate.
## when finished use scp or rsync to bring the tarballed .gz results file to your computer and open the .html file to evaluate the quality of your raw data.
#cd $CS
#tar cvzf $CS.gz $CS/*
#multiqc $CS/

######## FunGen Course Instructions ############
## Purpose: The purpose of this script is to trim sequencing adapters and low quality regions from the sequence read data with Trimmomatic,
##       and then use FASTQC to evaluate the quality of the data: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
## Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic
##              Input Data: Raw R1 & R2 reads (FASTQ); Adapter sequences to remove (FASTA)
##                              Downloaded read files, R1 and R2 files for each sample if paired-end data (FASTQ)
##              Output: Trimmed R1 & R2 paired and unpaired reads (FASTQ)
## FASTQC output is a folder for each file. The last line of this script will make a tarball of the output directory to bring back to your computer
##              Input Data: Raw R1 & R2 reads (FASTQ); Adapter sequences to remove (FASTA)
##                              Downloaded read files, R1 and R2 files for each sample if paired-end data (FASTQ)
##              Output: Trimmed R1 & R2 paired and unpaired reads (FASTQ)

## Purpose: The purpose of this script is to trim sequencing adapters and low quality regions from the read data.
## Input Data: Raw R1 & R2 reads (FASTQ); Adapter sequences to remove (FASTA)
## Output Data: Trimmed R1 & R2 paired and unpaired reads (FASTQ)
## More Information: http://www.usadellab.org/cms/?page=trimmomatic

#adapters=AdaptersToTrim_All.fa  ## This is a fasta file that has a list of adapters commonly used in NGS sequencing.
                                        ## You will likely need to edit this for other projects based on how your libraries
                                        ## were made to search for the correct adapters for your project

################ Trimmomatic ###################################
## Move to Raw Data Directory
#cd $DD

### Make list of file names to Trim
        ## this line is a set of piped (|) commands
        ## ls means make a list,
        ## grep means grab all the file names that end in ".fastq",
        ## cut that name into elements every where you see "_" and keep the first element (-f 1)
        ## sort the list and keep only the unique names and put it into a file named "list"
#ls | grep ".fq" |cut -d "_" -f 2 | sort | uniq > list


### Copy over the list of Sequencing Adapters that we want Trimmomatic to look for (along with its default adapters)
        ## CHECK: You may need to edit this path for the file that is in the class_shared directory from your account.
#cp /home/$MyID/class_shared/AdaptersToTrim_All.fa .
########这个地方可能需要修改， 目前的文件是从课上老师那里得到的

### Run a while loop to process through the names in the list and Trim them with the Trimmomatic Code
#while read i
#do

        ### Run Trimmomatic in paired end (PE) mode with 6 threads using phred 33 quality score format.
        ## STOP & DISCUSS: Check out the trimmomatic documentation to understand the parameters in line 77

#        java -jar /mnt/beegfs/home/aubmxa/.conda/envs/BioInfo_Tools/share/trimmomatic-0.39-1/trimmomatic.jar  PE -threads 6 -phred33 \
#        *"$i"_1.fq.gz *"$i"_2.fq.gz  \
#        $CD/"$i"_1_paired.fastq $CD/"$i"_1_unpaired.fastq  $CD/"$i"_2_paired.fastq $CD/"$i"_2_unpaired.fastq \
#        ILLUMINACLIP:AdaptersToTrim_All.fa:2:35:10 HEADCROP:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:6:30 MINLEN:36

                ## Trim read for quality when quality drops below Q30 and remove sequences shorter than 36 bp
                ## PE for paired end phred-score-type  R1-Infile   R2-Infile  R1-Paired-outfile R1-unpaired-outfile R-Paired-outfile R2-unpaired-outfile  Trimming paramenter
                ## MINLEN:<length> #length: Specifies the minimum length of reads to be kept.
                ## SLIDINGWINDOW:<windowSize>:<requiredQuality>  #windowSize: specifies the number of bases to average across
                ## requiredQuality: specifies the average quality required.
				
		#use trim_galore 
		trim_galore -output_dir $CD2 --phred -stringency 3 --paired $CD/"$i".fq.gz
        ############## FASTQC to assess quality of the Cleaned sequence data
        ## FastQC: run on each of the data files that have 'All' to check the quality of the data
        ## The output from this analysis is a folder of results and a zipped file of results

#fastqc $CD/"$i"_1_paired.fastq --outdir=$CS2
#fastqc $CD/"$i"_2_paired.fastq --outdir=$CS2

#fastqc $CD2/*_1 --outdir=$CS22
#fastqc $CD2/*_2 --outdir=$CS22

#done<list                       # This is the end of the loop
				
				
#########################  Now compress your results files from the Quality Assessment by FastQC
## move to the directory with the cleaned data
#cd $CS2

#######  Tarball the directory containing the FASTQC results so we can easily bring it back to our computer to evaluate.
## when finished use scp or rsync to bring the .gz file to your computer and open the .html file to evaluate
#tar cvzf $CS2.tar.gz $CS2/*				
#multiqc $CS2

#cd $CS22
#tar cvzf $CS22.tar.gz $CS22/*
#multiqc $CS22
				
######### FunGen Course Instructions ############
## Purpose: The purpose of this script is to 
##    Use HiSat2 to index your reference genome and then map your cleaned (paired) reads to the indexed reference
##              First need to use gffread to convert annotation file from .gff3 to .gft formate
##              Use Stringtie to count the reads mapped to genes and transcripts, defined in this case by the genome annotation file
##              use the python script to take the Stringtie results to make two counts matricies, one at the gene level and one at the transcript level
## HiSat2  Indexing  InPut: Reference genome file (.fasta), and annotation file (.gff3) (Optional)
##                    Output: Indexed genome 
## HiSat2 Mapping     Input: Cleaned read files, paired (.fasq); Indexed genome
##                    Output: Alignment .sam files  
## Samtools  Convert .sam to .bam and sort         Input: Alignment files  .sam
##                                                  Output: Sorted .bam files
## Stringtie  Counting reads  Input: sorted .bam file
##                            Output:  Directories of counts files for Ballgown (R program for DGE)
##              prepDE.py    Python script to create a counts matrics from the Stringtie output.  Inputs: Directory from Stringtie
##                                                                                                Output:  .csv files of counts matrix				
				
#  Set the stack size to unlimited
#ulimit -s unlimited

# Turn echo on so all commands are echoed in the output log
#set -x

 #This is what the "easy name" will be for the genome reference 



##################  Prepare the Reference Index for mapping with HiSat2   #############################
#cd $REFD_other
#cp ../../$REF_other.fasta .
#cp ../../$REF_other.gff .
#cp ../../$REF_other.gtf .
###  Identify exons and splice sites
#gffread $REF.gff3 -T -o $REF.gtf               ## gffread converts the annotation file from .gff3 to .gft formate for HiSat2 to use.
#extract_splice_sites.py $REF_other.gtf > $REF_other.ss
#extract_exons.py $REF_other.gtf > $REF_other.exon

#### Create a HISAT2 index for the reference genome. NOTE every mapping program will need to build a its own index.
#hisat2-build --ss $REF_other.ss --exon $REF_other.exon $REF_other.fasta ${REF_other}_index
#hisat2-build -p 10 $REF_other.fasta ${REF_other}_nonannotation_index

cd $REFD_origin
cp ../../$REF_origin.fna .
cp ../../$REF_origin.gff .

hisat2-build -p 10 $REF_origin.fna ${REF_origin}_index


				
# Move to the data directory
#cd $CD  #### This is where our clean paired reads are located.

## Create list of fastq files to map.    Example file format of your cleaned reads file names: SRR629651_1_paired.fastq SRR629651_2_paired.fastq
## grab all fastq files, cut on the underscore, use only the first of the cuts, sort, use unique put in list
#ls | grep ".fastq" |cut -d "_" -f 1| sort | uniq > list    #should list Example: SRR629651

## Move to the directory for mapping
#cd $MAPD_other

## move the list of unique ids from the original files to map
#cp $CD/list  . 

#while read i;
#do
  ## HiSat2 is the mapping program
  ##  -p indicates number of processors, --dta reports alignments for StringTie --rf is the read orientation
#   hisat2 -p 6 --dta --phred33       \
#    -x "$REFD_other"/${REF_other}_index       \
#    -1 "$CD"/"$i"_1_paired.fastq  -2 "$CD"/"$i"_2_paired.fastq      \
#    -S "$i".sam

    ### view: convert the SAM file into a BAM file  -bS: BAM is the binary format corresponding to the SAM text format.
    ### sort: convert the BAM file to a sorted BAM file.
    ### Example Input: SRR629651.sam; Output: SRR629651_sorted.bam
#  samtools view -@ 6 -bS "$i".sam > "$i".bam  ### This works on ASC

    ###  This is sorting the bam
#  samtools sort -@ 6  "$i".bam    "$i"_sorted

    ### Index the BAM and get mapping statistics
#  samtools flagstat   "$i"_sorted.bam   > "$i"_Stats.txt

  ### Stringtie is the program that counts the reads that are mapped to each gene, exon, transcript model. 
  ### The output from StringTie are counts folders in a directory that is ready to bring into the R program Ballgown to 
  ### Original: This will make transcripts using the reference geneome as a guide for each sorted.bam
  ### eAB options: This will run stringtie once and  ONLY use the Ref annotation for counting readsto genes and exons 
  ###
#mkdir "$COUNTSD_other"/"$i"
#stringtie -p 6 -e -B -G  "$REFD_other"/"$REF_other".gff -o "$COUNTSD_other"/"$i"/"$i".gtf -l "$i"   "$MAPD_other"/"$i"_sorted.bam

#done<list

#####################  Copy Results to home Directory.  These will be the files you want to bring back to your computer.
### these are your stats files from Samtools
#cp *.txt $RESULT_other

### The prepDE.py is a python script that converts the files in your ballgown folder to a count matrix
#cd $COUNTSD_other
#python /home/$MyID/class_shared/prepDE.py -i $COUNTSD_other

### copy the final results files (the count matricies that are .cvs to your home directory)
#cp *.csv $RESULT_other







## Move to the directory for mapping
#cd $MAPD_other_nonannotation

## move the list of unique ids from the original files to map
#cp $CD/list  . 

#while read i;
#do
  ## HiSat2 is the mapping program
  ##  -p indicates number of processors, --dta reports alignments for StringTie --rf is the read orientation
#   hisat2 -p 6 --dta --phred33       \
#    -x "$REFD_other"/${REF_other}_nonannotation_index       \
#    -1 "$CD"/"$i"_1_paired.fastq  -2 "$CD"/"$i"_2_paired.fastq      \
#    -S "$i".sam

    ### view: convert the SAM file into a BAM file  -bS: BAM is the binary format corresponding to the SAM text format.
    ### sort: convert the BAM file to a sorted BAM file.
    ### Example Input: SRR629651.sam; Output: SRR629651_sorted.bam
#  samtools view -@ 6 -bS "$i".sam > "$i".bam  ### This works on ASC

    ###  This is sorting the bam
#  samtools sort -@ 6  "$i".bam    "$i"_sorted

    ### Index the BAM and get mapping statistics
#  samtools flagstat   "$i"_sorted.bam   > "$i"_Stats.txt

  ### Stringtie is the program that counts the reads that are mapped to each gene, exon, transcript model. 
  ### The output from StringTie are counts folders in a directory that is ready to bring into the R program Ballgown to 
  ### Original: This will make transcripts using the reference geneome as a guide for each sorted.bam
  ### eAB options: This will run stringtie once and  ONLY use the Ref annotation for counting readsto genes and exons 
  ###
#mkdir "$COUNTSD_other_nonannotation"/"$i"
#stringtie -p 6 -e -B -G  "$REFD_other"/"$REF_other".gff -o "$COUNTSD_other_nonannotation"/"$i"/"$i".gtf -l "$i"   "$MAPD_other_nonannotation"/"$i"_sorted.bam

#done<list

#####################  Copy Results to home Directory.  These will be the files you want to bring back to your computer.
### these are your stats files from Samtools
#cp *.txt $RESULT_other_nonannotation

### The prepDE.py is a python script that converts the files in your ballgown folder to a count matrix
#cd $COUNTSD_other_nonannotation
#python /home/$MyID/class_shared/prepDE.py -i $COUNTSD_other_nonannotation

### copy the final results files (the count matricies that are .cvs to your home directory)
#cp *.csv $RESULT_other_nonannotation





## Move to the directory for mapping
cd $MAPD_origin

## move the list of unique ids from the original files to map
cp $CD/list  . 

while read i;
do
  ## HiSat2 is the mapping program
  ##  -p indicates number of processors, --dta reports alignments for StringTie --rf is the read orientation
   hisat2 -p 6 --dta --phred33       \
    -x "$REFD_origin"/${REF_origin}_index       \
    -1 "$CD"/"$i"_1_paired.fastq  -2 "$CD"/"$i"_2_paired.fastq      \
    -S "$i".sam

    ### view: convert the SAM file into a BAM file  -bS: BAM is the binary format corresponding to the SAM text format.
    ### sort: convert the BAM file to a sorted BAM file.
    ### Example Input: SRR629651.sam; Output: SRR629651_sorted.bam
  samtools view -@ 6 -bS "$i".sam > "$i".bam  ### This works on ASC

    ###  This is sorting the bam
  samtools sort -@ 6  "$i".bam    "$i"_sorted

    ### Index the BAM and get mapping statistics
  samtools flagstat   "$i"_sorted.bam   > "$i"_Stats.txt

  ### Stringtie is the program that counts the reads that are mapped to each gene, exon, transcript model. 
  ### The output from StringTie are counts folders in a directory that is ready to bring into the R program Ballgown to 
  ### Original: This will make transcripts using the reference geneome as a guide for each sorted.bam
  ### eAB options: This will run stringtie once and  ONLY use the Ref annotation for counting readsto genes and exons 
  ###
mkdir "$COUNTSD_origin"/"$i"
stringtie -p 6 -e -B -G  ../../lizard.gff -o "$COUNTSD_origin"/"$i"/"$i".gtf -l "$i"   "$MAPD_origin"/"$i"_sorted.bam

done<list

#####################  Copy Results to home Directory.  These will be the files you want to bring back to your computer.
### these are your stats files from Samtools
cp *.txt $RESULT_origin

### The prepDE.py is a python script that converts the files in your ballgown folder to a count matrix
cd $COUNTSD_origin
python /home/$MyID/class_shared/prepDE.py -i $COUNTSD_origin

### copy the final results files (the count matricies that are .cvs to your home directory)
cp *.csv $RESULT_origin
