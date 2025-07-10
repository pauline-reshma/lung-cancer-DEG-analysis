# Load required libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(apeglm)

# STEP 1: Load data
count_data <- read.csv("../data/counts_matrix_deseq2_ready.csv", row.names = 1)
col_data <- read.csv("../data/sample_metadata_deseq2_ready.csv", row.names = 1)

# STEP 2: Clean NA values
count_data[is.na(count_data)] <- 0

# STEP 3: Prepare DESeq2 dataset
col_data$Condition <- as.factor(col_data$Condition)
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = ~ Condition)

# STEP 4: Filter low-expression genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# STEP 5: Run DESeq2
dds <- DESeq(dds)

# STEP 6: Get results
res <- results(dds, alpha = 0.05)

# STEP 7: Shrink log2 Fold Change
res <- lfcShrink(dds, coef = 2, res = res)

# STEP 8: Save DEGs (filtered)
resOrdered <- res[order(res$padj), ]
res_filtered <- subset(resOrdered, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(res_filtered), "../results/DEGs_filtered_FDR0.05_Log2FC1.csv")

# STEP 9: Volcano plot
res$threshold <- as.factor(abs(res$log2FoldChange) > 1 & res$padj < 0.05)
ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Expression", x = "Log2 Fold Change", y = "-log10(FDR)") +
  ggsave("../results/volcano_plot.png")

# STEP 10: Heatmap (Top 50 DEGs)
select_genes <- head(order(res$padj), 50)
vsd <- vst(dds, blind = FALSE)
pheatmap(assay(vsd)[select_genes, ], cluster_rows = TRUE, show_rownames = TRUE,
         cluster_cols = TRUE, annotation_col = col_data)
