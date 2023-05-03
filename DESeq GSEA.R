# Load DESeq2 results file
results <- read.csv("Heart_Liver_DESeq2.csv", stringsAsFactors = FALSE)
#results <- read.csv("Heart_Muscle_DESeq2.csv", stringsAsFactors = FALSE)
#results <- read.csv("Liver_Muscle_DESeq2.csv", stringsAsFactors = FALSE)

# Remove rows with "-Inf" log2fold change values
#results <- results[!is.infinite(results$log2FoldChange) | results$log2FoldChange > 0,]

# Replace "-Inf" values with very small or very large values
#results$log2FoldChange[results$log2FoldChange == -Inf] <- -1000
#results$log2FoldChange[results$log2FoldChange == Inf] <- 1000

# Extract gene names, log2fold change, and p-value
genes <- results$gene_id
logFC <- results$log2FoldChange
pvalue <- results$pvalue

# Replace any "-Inf" values in logFC with zero
#logFC[is.infinite(logFC) & logFC < 0] <- 0

# Calculate the rank based on log2fold change and p-value
rank <- ifelse(logFC > 0, -log10(pvalue), log10(pvalue))
names(rank) <- genes
rank <- rank[order(rank, decreasing = TRUE)]

# Save the ranked list as a tab-delimited file with the .rnk extension
write.table(rank, file = "DESeq2_ranked_list_Heart_Liver.rnk", sep = "\t", col.names = FALSE)
#write.table(rank, file = "DESeq2_ranked_list_Heart_Muscle.rnk", sep = "\t", col.names = FALSE)
#write.table(rank, file = "DESeq2_ranked_list_Liver_Muscle.rnk", sep = "\t", col.names = FALSE)
