#!/bin/bash

# Extracting raw reads count. The following script will read through the individual column files & 
# Will grep the 1st and 4th column of the count files and will generate new files with gene names & read numbers.


## Create a directory that will contain the outputs
mkdir ./full_data/deseq2/counts_4thcol


for i in ./full_data/counts_star/*ReadsPerGene.out.tab
do echo $i
# retrieve the first (gene name) and fourth column (raw reads)
cut -f1,4 $i | grep -v "_" > ./full_data/deseq2/counts_4thcol/`basename $i ReadsPerGene.out.tab`_counts.txt
done

