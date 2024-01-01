## Differential gene expression and Gene Set Enrichment Analysis
After Step 6 the analysis is continued on R-programming language. We will use [DeSeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and [fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html) packages.

**Set the working directory to your current project directory**
```
setwd("YOUR/WORKING/DIRECTORY/PATH")
```
**Load necessary libraries** 
```
library(DESeq2)
library(tidyverse)
library(apeglm)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(msigdbr)
library(magrittr)
library(fgsea)
library(ggridges)
library(enrichplot)
library(DOSE)
```
**Reads in count and metadata from user-specified files**
 - countData : count dataframe
 - colData : sample metadata in the dataframe with row names as sampleID's
 - design : The design of the comparisons to use. 
           Use (~) before the name of the column variable to compare
```
counts_data <- read.table(counts_data_file, header = TRUE, row.names = 1)
colData <- read.table(metadata_file, header = TRUE, row.names = 1)
```
Check if all column names in counts_data match row names in colData
```
stopifnot(all(colnames(counts_data) %in% rownames(colData)))
stopifnot(all(colnames(counts_data) == rownames(colData)))

```
**Create DESeq2 dataset**

```
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData, 
                              design = ~ condition)
```
**Filter low count genes (my sample size was 10)**  
```
dds <- dds[rowSums(counts(dds)) >= 9,] 
```
**Run statistical analysis**
```
dds <- DESeq(dds)
#Get results from testing with FDR adjust pvalues
res <- results(dds, pAdjustMethod = "fdr", alpha = 0.05))
```

**Add gene symbols to results**
```
res$symbol <- mapIds(org.Mm.eg.db, 
                     keys = gsub("\\..*","",rownames(res)), 
                     column = "SYMBOL", 
                     keytype = "ENSEMBL", 
                     multiVals = "first")
# Add gene symbol
res$symbol <- row.names(res)

# Add ENTREZ ID
res$entrez <- mapIds(x = org.Mm.eg.db,
                     keys = row.names(res),
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")

# Subset for only significant genes (q < 0.05)
res_sig <- subset(res, padj < 0.05)
head(res_sig)
```
**Write results to .txt files**
```
# Write normalized gene counts to a .txt file
write.table(x = as.data.frame(counts(dds), normalized = T), 
            file = 'DESeq2_normalized_counts.txt', 
            sep = '\t', 
            quote = F,
            col.names = NA)

# Write significant normalized gene counts to a .txt file
write.table(x = counts(dds[row.names(res_sig)], normalized = T), 
            file = 'DESeq2_normalized_counts_significant.txt', 
            sep = '\t', 
            quote = F, 
            col.names = NA)

# Write the annotated results table to a .txt file
write.table(x = as.data.frame(res), 
            file = "DESEq2_results_gene_annotated.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)

# Write significant annotated results table to a .txt file
write.table(x = as.data.frame(res_sig), 
            file = "DESEq2_results_gene_annotated_significant.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)
```

**Plot the results**
- MA plot 
```
plotMA(res, ylim=c(-2,2))
```
- PCA plot 
```
# Transform data using regularized log (rlog)
rld <- rlog(dds, blind = FALSE)

# Perform PCA and store the results
pcaData <- plotPCA(rld, intgroup = "condition", returnData = TRUE)

# Calculate the percentage of variance explained by PC1 and PC2
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Generate PCA plot using ggplot2
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape = condition)) +
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  geom_text(aes(label = rownames(pcaData)), 
            nudge_x = 0.25, nudge_y = 1.00, 
            check_overlap = TRUE) +
  coord_fixed()
``` 
- Heatmap for Top Genes
```
# Load your data 
num_top_genes <- 25  # Number of top genes to display

# Select top genes based on adjusted p-value and sort by log2FoldChange
top_genes <- dge_data %>%
             dplyr::arrange(desc(log2FoldChange)) %>%
             head(num_top_genes)

# Prepare matrix for heatmap
mat <- assay(rld)[rownames(top_genes), colnames(dge_data)]
mat_scaled <- t(apply(mat, 1, scale))  # Z-score scaling

# Create and plot heatmap
Heatmap(mat_scaled, name = "Z-score", 
        cluster_rows = FALSE, 
        cluster_columns = TRUE)
```
- Volcano Plot 
```
# Define thresholds and load data
log2FC_threshold_up <- 2
log2FC_threshold_down <- -2
pvalue_threshold <- 0.05
dge_data_file <- "path/to/your/dge_data.csv"  # Replace with your file path
dge_data <- read.csv(dge_data_file)

# Categorize genes and prepare for plotting
dge_data$Category <- ifelse(dge_data$log2FoldChange > log2FC_threshold_up & dge_data$padj < pvalue_threshold, 'Upregulated',
                           ifelse(dge_data$log2FoldChange < log2FC_threshold_down & dge_data$padj < pvalue_threshold, 'Downregulated', 'Not Significant'))
labeled_genes <- subset(dge_data, Category %in% c('Upregulated', 'Downregulated'))

# Create the volcano plot
volcano_plot <- ggplot(dge_data, aes(x = log2FoldChange, y = -log10(padj), color = Category)) +
                geom_point(alpha = 0.5) +
                geom_vline(xintercept = c(log2FC_threshold_up, log2FC_threshold_down), linetype = "dashed", color = "black", alpha = 0.5) +
                geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "black", alpha = 0.5) +
                scale_color_manual(values = c('Not Significant' = 'grey50', 'Upregulated' = 'red', 'Downregulated' = 'blue')) +
                theme_minimal() +
                labs(title = "Volcano Plot: Differential Gene Expression", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
                coord_cartesian(xlim = c(-10, 10), ylim = c(0, -log10(pvalue_threshold) + 5)) +
                geom_text_repel(data = labeled_genes, aes(label = Gene), max.overlaps = Inf)

# Print the plot
print(volcano_plot)

```
- Visualization of Gene of Interest
 ```
gene_id <- "ENS....." # Replace with your specific gene ID
# Prepare expression data for plotting
tcounts <- as.data.frame(log2(counts(dds, normalized = TRUE)[gene_id, ] + 0.5))
tcounts$condition <- colData(dds)$condition

# Create the plot using ggplot2
ggplot(tcounts, aes(x = condition, y = gene_id, fill = condition)) + 
  geom_point() +
  ylim(0, 18) +
  labs(x = "Condition",
       y = "Expression (log normalized counts)",
       title = paste("Expression of gene", gene_id))
```

**GeneSetEnrichmentAnalysis (KEGG) and Visualisation**

The same method can also be applied using gseGO for getting results based on the Gene Ontology database.

```
kk2 <- gseKEGG(geneList = gene_list, 
               organism = "mmu", 
               nPerm = 10000,
               minGSSize = 3,
               maxGSSize = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               keyType = "ncbi-geneid")

# Dotplot visualization of top 10 enriched categories
dotplot_gsea <- dotplot(kk2, showCategory = 10) + facet_grid(. ~ .sign)

# Enhanced matrix plot showing relationships between terms
emapplot_gsea <- emapplot(pairwise_termsim(kk2), showCategory = 20)

# Convert enrichment results to a more readable format
edox <- setReadable(kk2, 'org.Mm.eg.db', 'ENTREZID')

# Create network plots with different node labels
p1 <- cnetplot(edox, node_label = "category", cex_label_category = 1.2, foldChange = geneList)
p2 <- cnetplot(edox, node_label = "gene", cex_label_gene = 0.8, foldChange = geneList)
p3 <- cnetplot(edox, node_label = "all", foldChange = geneList)

# Combine plots into a single grid for comparison
combined_cnet_plots <- plot_grid(p1, p2, p3, ncol = 3, labels = "AUTO")

# Print the combined network plots
print(combined_cnet_plots)
```
