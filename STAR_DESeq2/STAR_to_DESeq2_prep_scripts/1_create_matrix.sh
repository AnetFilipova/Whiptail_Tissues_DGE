#!/bin/bash

# Create a matrix with the counts
cd ./full_data/counts_star

# retrieve the 4th column of each "ReadsPerGene.out.tab" file + the first column that contains the gene IDs
paste *ReadsPerGene.out.tab | grep -v "_" | awk '{printf "%s\t", $1}{for (i=4;i<=NF;i+=4) printf "%s\t", $i; printf "\n" }' > tmp

# add header: "gene_name" + the name of each of the counts file
sed -e "1igene_name\t$(ls *ReadsPerGene.out.tab | tr '\n' '\t' | sed 's/ReadsPerGene.out.tab//g')" tmp | cut -f1-7 > ../deseq2/raw_counts_matrix.txt

# another way can be the following one
ls *.tab | awk 'BEGIN{ORS="";print "gene name\t"}{print $0"\t"}END{print "\n"}'| sed 's/ReadsPerGene.out.tab//g' > ../deseq2/raw_counts_matrix.txt; cat tmp >> ../deseq2/raw_counts_matrix.txt


# remove temporary file
rm tmp