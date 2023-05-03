# Load edgeR results file
results <- read.csv("Heart_vs_Liver_DEG.csv", stringsAsFactors = FALSE)
#results <- read.csv("Heart_vs_Muscle_DEG.csv", stringsAsFactors = FALSE)
#results <- read.csv("Liver_vs_Muscle_DEG.csv", stringsAsFactors = FALSE)

# Extract gene names, log2fold change, and p-value
genes <- results$GeneID
logFC <- results$logFC
pvalue <- results$PValue

# Calculate the rank based on log2fold change and p-value
rank <- ifelse(logFC > 0, -log10(pvalue), log10(pvalue))
names(rank) <- genes
rank <- rank[order(rank, decreasing = TRUE)]

# Save the ranked list as a tab-delimited file with the .rnk extension
write.table(rank, file = "edgeR_ranked_list_Heart_Liver.rnk", sep = "\t", col.names = FALSE)
#write.table(rank, file = "edgeR_ranked_list_Heart_Muscle.rnk", sep = "\t", col.names = FALSE)
#write.table(rank, file = "edgeR_ranked_list_Liver_Muscle.rnk", sep = "\t", col.names = FALSE)

