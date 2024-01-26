# Loading the library for DEG
library(DESeq2)
# Loading the library for plots
library("ggplot2")
# Loading the library for volcano plot
library("EnhancedVolcano")
# Loading the library for heatmap
library("pheatmap")

# # Setting the Path to the featureCounts Count Matrix
# setwd("RNA-seq-protocol/quantification_featureCounts/")
pdf(file="PlotsDifferentialExpression.pdf")

# Reading the Count Matrix
countData <- read.delim("matriz_arabidopsis_2023_gen.txt", header =TRUE, row.names = 1, skip=1)

extracted_names <- sub(".+\\.(\\w+)\\..+$", "\\1", colnames(countData))
# colnames(your_data) <- new_column_names

sample_dictionary <- c("SRR10207204" = "Control", "SRR10207210" = "Control", "SRR10207216" = "Control", "SRR10207206" = "Treatment", "SRR10207212" = "Treatment", "SRR10207218" = "Treatment")
# Add more entries as needed

# Map the sample types using the dictionary
sample_types <- sapply(extracted_names, function(name) sample_dictionary[[name]])

# Combine the extracted names with the sample types
# new_column_names <- paste(extracted_names, sample_types, sep = "_")

# Replace the column names
# colnames(your_data) <- new_column_names

# Replace the column names
colnames(countData) <- extracted_names

# Description of Samples
# condition <- factor(c("Control", "Treatment", "Control", "Treatment", "Control", "Treatment"))
# colData <- data.frame(row.names = colnames(countData), condition)
# sample_types <- sapply(extracted_names, function(name) sample_dictionary[[name]])

# Create the condition factor using the sample types
condition <- factor(sample_types)
colData <- data.frame(row.names = colnames(countData), condition)

# print(head(colData))

# Creating a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~condition)

# Generating a PCA Plot
rld <- rlog(dds, blind = F)
# pca_plot <- plotPCA(rld, intgroup = "condition") + geom_text(aes(label=name),vjust=0.2)
plotPCA(rld, intgroup = "condition") + geom_text(aes(label=name),vjust=0.2)
# ggsave("PCA_plot.png", pca_plot, dpi = 400)
# ggsave("PCA_plot.png", pca_plot, dpi = 400)

# Filtering Genes with Very Low Expression (less than 10 reads)
dds <- dds[rowSums(counts(dds)) > 10, ]

# Performing Differential Expression Analysis
dds <- DESeq(dds)

# Extracting Results
res <- results(dds)

# Sorting the Summary List by Adjusted p-adj
res <- res[order(res$padj), ]

# Extracting Contrasts between Conditions
Treatment_vs_Control <- results(dds, contrast = c("condition","Treatment", "Control"))
# write.csv(Treatment_vs_Control, file = "all.csv")
# print(Treatment_vs_Control)

# Obtaining a List of Differentially Expressed Genes with Stricter Filtering
deg <- subset(Treatment_vs_Control, padj < 0.05 & abs(log2FoldChange) > 1)
# print(deg)
# Exporting the DEG Table
write.csv(deg, file = "DEG_Treatment_vs_Control_strict.csv")
# Printing and Exporting up-DEGs
up <- subset(deg, log2FoldChange > 1)
write.csv(up, file = "Up_Treatment_vs_Control_strict.csv")
# Printing and Exporting down-DEGs
down <- subset(deg, log2FoldChange < (-1))
write.csv(down, file = "Down_Treatment_vs_Control_strict.csv")

# Generating a Volcano Plot with EnhancedVolcano
EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'pvalue')
# volcano_plot <- EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'pvalue')
# ggsave("volcano_plot.png", volcano_plot, dpi = 400)

# Setting the Contrast of Interest
# plotma <- plotMA(Treatment_vs_Control, alpha = 0.05, main = "Inoculated with Trichoderma vs Control", xlab = "mean of normalized counts",  ylim = c(-2, 2))
plotMA(Treatment_vs_Control, alpha = 0.05, main = "Inoculated with Trichoderma vs Control", xlab = "mean of normalized counts",  ylim = c(-2, 2))
# plotma <- plotMA(Treatment_vs_Control, alpha = 0.5, main = "Inoculated with Trichoderma vs Control", xlab = "mean of normalized counts")
# plotma <- plotMA(res, alpha = 0.05, main = "Inoculated with Trichoderma vs Control", xlab = "mean of normalized counts")
# plotma <- plotMA(Treatment_vs_Control, ylim = c(-2, 2))
# ggsave("MA_plot.png", plotma, dpi = 400)

# Selecting the Top 20 Differentially Expressed Genes in the Treatment_vs_Control Comparison
res_ordered <- Treatment_vs_Control[order(Treatment_vs_Control$padj), ]
# top_genes <- row.names(res_ordered)[1:20]
top_genes <- row.names(res_ordered)[1:10]

# Extracting and Normalizing Counts
counts <- counts(dds, normalized = TRUE)
counts_top <- counts[top_genes, ]

# Applying a Logarithmic Transformation to Counts
log_counts_top <- log2(counts_top + 1)

# Creating an Annotation Data Frame Based on Condition Information (colData)
df <- colData

# Generating a Heatmap using the pheatmap library
pheatmap(log_counts_top, annotation = df)
# heatmap_20 <- pheatmap(log_counts_top, annotation = df)
# ggsave("heatmap_plot.png", heatmap_20, dpi = 400)
