#!/bin/sh



#! /bin/bash

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
# Loading modules
module load star/
#module load samtools/1.11
module load stringtie/2.2.1
module load gffcompare
module load python/2.7.1
module load gcc/9.3.0
module load samtools
module load bcftools/1.2
module load gffread/
module load gffcompare/


#  Set the stack size to unlimited
ulimit -s unlimited

# Turn echo on so all commands are echoed in the output log
set -x



MyID=aubclsb0312

WD=/scratch/$MyID/Whiptail_Samples  
CLEAND=/scratch/aubclsb0312/Whiptail_Samples/CleanData
REFD=/scratch/aubclsb0312/Whiptail_Samples/ref/PodMur                                              ## Example:/scratch/$MyID/PracticeRNAseq/DaphniaRefGenome_6                      # this directory contains the indexed reference genome for the garter snake
MAPD=/scratch/$MyID/Whiptail_Samples/Map_STAR                                             
COUNTSD=/scratch/$MyID/Whiptail_Samples/Counts_StringTie_6

RESULTSD=/scratch/$MyID/Whiptail_Samples/Counts_STAR_ST_6      

REF=GCF_004329235.1_PodMur_1.0_genomic

#mkdir -p $REFD
mkdir -p $MAPD
mkdir -p $COUNTSD
mkdir -p $RESULTSD


cd $REFD



### Generating genome indices
STAR --runThreadN 24
--runMode genomeGenerate
--genomeDir ${REFD}
--genomeFastaFiles ${REFD}/${REF}.fna
--sjdbGTFfile ${REFD}/${REF}.gtf
--sjdbOverhang ReadLength-1



cd $CLEAND

ls | grep ".fastq.gz" |cut -d "_" -f 1| sort | uniq > list


#cd $MAPD

#mv $CLEAND/list  .


while read i;
do


STAR  --genomeDir ${REFD} \
--readFilesIn "$i"_1_paired.fastq.gz "$i"_2_paired.fastq.gz \
--readFilesCommand zcat --runThreadN 24 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix alignments_STAR/"$i"



done<list
