# The script can be used to statistically analyze differential expression of genes
# and generate visualizations of the DEGs using outputs of Star mapped files
# The inputs used here were pre-processed using bash scripts found in the main directory

# The following script was modified from https://biocorecrg.github.io/RNAseq_course_2019/differential_expression.html


# Load the DESeq2 package

library(DESeq2)

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
nrow(sampletable) # if this is not 6, please raise your hand !
ncol(sampletable) # if this is not 4, also raise your hand !

# create the DESeq object
    # countData is the matrix containing the counts
    # sampletable is the sample sheet / metadata we created
    # design is how we wish to model the data: what we want to measure here is the difference between the treatment times
# Option 1 that reads in a matrix (we will not do it here):

# first read in the matrix
count_matrix <- read.delim("raw_counts_matrix.txt", header=T, sep="\t", row.names=1)
colnames(count_matrix) <- sampletable$SampleName # The sample names are not showing in columns
head(count_matrix)
ncol(count_matrix)
nrow(count_matrix)

	# then create the DESeq object
		# countData is the matrix containing the counts
		# sampletable is the sample sheet / metadata we created
		# design is how we wish to model the data: what we want to measure here is the difference between the treatment times


se_star <- DESeqDataSetFromMatrix(countData = count_matrix,
                                  colData = sampletable,
                                  design = ~TISSUE)


# -----------------  F I L T E R I N G -------------------------
# Number of genes before filtering:
nrow(se_star)

# Filter
se_star <- se_star[rowSums(counts(se_star)) > 10, ]


# Number of genes left after low-count filtering:
nrow(se_star)

# ##################  Fitting Statistical Model #################################

se_star2 <- DESeq(se_star)
head(se_star2$SampleName)

# ##################  Extracting Results #################################

# compute normalized counts (log2 transformed); + 1 is a count added to avoid errors during the log2 transformation: log2(0) gives an infinite number, but log2(1) is 0.
# normalized = TRUE: divide the counts by the size factors calculated by the DESeq function
norm_counts <- log2(counts(se_star2, normalized = TRUE)+1)
head(norm_counts)

# read in the file containing gene symbols
tx2gene <- read.table("annotation.csv", 
                      sep="\t",
                      header=F)
head(tx2gene)

# add the gene symbols
norm_counts_symbols <- merge(unique(tx2gene[,2:3]), data.frame(ID=rownames(norm_counts), norm_counts), by=1, all=F)
head(norm_counts_symbols)

# To bring the lost samples names back
colnames(norm_counts) <- se_star2$SampleName
head(norm_counts)


# write normalized counts to text file
write.table(norm_counts_symbols, "normalized_counts.txt", quote=F, col.names=T, row.names=F, sep="\t")

# visualization

# # load libraries pheatmap to create the heatmap plot
library(pheatmap)

# Try with the vst transformation
vst <- vst(se_star2)

# calculate between-sample distance matrix
sampleDistMatrix <- as.matrix(dist(t(assay(vst))))

# create figure in PNG format
pheatmap(sampleDistMatrix, annotation_col = sampletable)

# png("PCA_star.png")
plotPCA(object = vst,
        intgroup = "TISSUE")

## check results names: depends on what was modeled. Here it was the "TISSUE"
resultsNames(se_star2)
se_star2$TISSUE

# extract results for t25 vs t0
# contrast: the column from the metadata that is used for the grouping of the samples (Time), then the baseline (t0) and the group compared to the baseline (t25) -> results will be as "t25 vs t0"
LiverHeart <- results(object = se_star2, 
              name="TISSUE_Liver_vs_Heart")
MuscleHeart <- results(object = se_star2, 
                       name="TISSUE_Liver_vs_Heart")

# processing the same results as above but including the log2FoldChange shrinkage
# useful for visualization and gene ranking
LH_shrink <- lfcShrink(dds = se_star2,
                       coef="TISSUE_Liver_vs_Heart",
                       type="apeglm")
MH_shrink <- lfcShrink(dds = se_star2,
                    coef="TISSUE_Liver_vs_Heart",
                    type="apeglm")

# check first rows of both results
head(LH_shrink)
head(MH_shrink)

# add the more comprehensive gene symbols to de_shrink
LH_symbols <- merge(unique(tx2gene[,2:3]), data.frame(ID=rownames(LH_shrink), LH_shrink), by=1, all=F)
head(LH_symbols)

MH_symbols <- merge(unique(tx2gene[,2:3]), data.frame(ID=rownames(LH_shrink), LH_shrink), by=1, all=F)
head(MH_symbols)

# write differential expression analysis result to a text file
write.table(LH_symbols, "liver_heart_deseq2_results.txt", quote=F, col.names=T, row.names=F, sep="\t")

write.table(MH_symbols, "muscle_heart_deseq2_results.txt", quote=F, col.names=T, row.names=F, sep="\t")


