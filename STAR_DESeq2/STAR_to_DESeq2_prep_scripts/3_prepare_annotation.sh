#!/bin/bash

cd ./full_data/deseq2

# Download annotation for all chromosomes
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz

# I downloaded the annotation of the latest release 39 (GRCh39) from the ENSEMBL website
# instead of the regular one (depends on the version of the genome used by Sheniqua for 
# the alignment using STAR mapper 

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/64176/100/GCF_004329235.1_PodMur_1.0/GCF_004329235.1_PodMur_1.0_genomic.gtf.gz

# first column is the transcript ID, second column is the gene ID, third column is the gene symbol
zcat GCF_004329235.1_PodMur_1.0_genomic.gtf.gz | awk -F "\t" 'BEGIN{OFS="\t"}{if($3=="transcript"){split($9, a, "\""); print a[4],a[2],a[8]}}' > tx2gene.gencode.v39.csv
