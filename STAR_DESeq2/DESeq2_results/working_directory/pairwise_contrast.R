# The script can be used to statistically analyze differential expression of genes
# and generate visualizations of the DEGs using outputs of Star mapped files
# The inputs used here were pre-processed using bash scripts found in the main directory

# The following script was modified from https://biocorecrg.github.io/RNAseq_course_2019/differential_expression.html


# Load the DESeq2 package

library(DESeq2)
library(tidyverse)
# read in the sample sheet
# header = TRUE: the first row is the "header", i.e. it contains the column names.
# sep = "\t": the columns/fields are separated with tabs.

# Copy the SampleSheet to the current directory
sampletable <- read.table("sample_sheet.txt", header=T, sep="\t")

# add the sample names as row names (it is needed for some of the DESeq functions)
rownames(sampletable) <- sampletable$SampleName

# display the first 6 rows
head(sampletable)

# check the number of rows and the number of columns
nrow(sampletable)
ncol(sampletable) 

# create the DESeq object
# countData is the matrix containing the counts
# sampletable is the sample sheet / metadata we created
# design is how we wish to model the data: what we want to measure here is the difference between the treatment times
# Option 1 that reads in a matrix (we will not do it here):

# first read in the matrix
count_matrix <- read.table("raw_counts_matrix.txt", header = T, row.names = 1)
head(count_matrix)
colnames(count_matrix) <- sampletable$SampleName # The sample names are not showing in columns
head(count_matrix)
ncol(count_matrix)
nrow(count_matrix)


# then create the DESeq object
# countData is the matrix containing the counts
# sampletable is the sample sheet / metadata we created
# design is how we wish to model the data: what we want to measure here is the difference between the treatment times

Deseq_contrast_data <- DESeqDataSetFromMatrix(countData = count_matrix,
                                  colData = sampletable,
                                  design = ~TISSUE)



# -----------------  F I L T E R I N G -------------------------
# Number of genes before filtering:
nrow(Deseq_contrast_data)

# Filter
Deseq_contrast_data <- Deseq_contrast_data[rowSums(counts(Deseq_contrast_data)) > 0, ]


# Number of genes left after low-count filtering:
nrow(Deseq_contrast_data)

# test contrast 

Deseq_con_MLE <- DESeq(Deseq_contrast_data, modelMatrixType="expanded", betaPrior=T)
Deseq_con_MLE


# Checking group names
resultsNames(Deseq_con_MLE)

# Create pairwise contrast matrix for all groups in the treatment column

contrast1 = c("TISSUE", "Liver", "Skeletal_Muscle")
contrast2 = c("TISSUE", "Heart", "Skeletal_Muscle")
contrast3 = c("TISSUE", "Heart", "Liver")

# running contrasts
contrast1_res_MLE <- results(Deseq_con_MLE, contrast=contrast1, alpha = 0.05)
contrast2_res_MLE <- results(Deseq_con_MLE, contrast=contrast2, alpha = 0.05)
contrast3_res_MLE <- results(Deseq_con_MLE, contrast=contrast3, alpha = 0.05)


# writing csv of all results
write.csv(as.matrix(contrast1_res_MLE), file = "Liver_Muscle_DESeq2.csv", row.names = T)
write.csv(as.matrix(contrast2_res_MLE), file = "Heart_Muscle_DESeq2.csv", row.names = T)
write.csv(as.matrix(contrast3_res_MLE), file = "Heart_Liver_DESeq2.csv", row.names = T)







