# RNA-seq-protocol
 This is the protocol for conducting RNAseq analysis in plants. This repository has been created to make it easier to follow the protocol presented in the Current Protocols journal.
## Protocol 1

### Raw Data for the analysis:

In the 7th step of the protocol outlined in Current Protocols, we provide detailed instructions for conducting a differential expression analysis. To reach this stage, it is imperative to download the Fastq files and complete the entire analysis, culminating in the contingency table containing counts for each gene. However, to simplify the initiation of hands-on experience with statistical analyses in R, we have deposited the final count matrix here for your convenience.

Matrix: [Matrix_arabidopsis_2023](https://raw.githubusercontent.com/jmvillalobos/RNA-seq-protocol/main/matriz_arabidopsis_2023.txt)


#Supplementary 2
#Differential gene expression analysis

#a.Library Installation
# Installing DESeq2
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

# Installing ggplot2
install.packages("ggplot2")

# Installing EnhancedVolcano
BiocManager::install("EnhancedVolcano")

# Installing pheatmap
install.packages("pheatmap")


#b. Loading the library for DEG
library(DESeq2)
# Loading the library for plots
library("ggplot2")
# Loading the library for volcano plot
library("EnhancedVolcano")
# Loading the library for heatmap
library("pheatmap")

#c.Setting the Path to the featureCounts Count Matrix
setwd("C:/project_2023/quantification_featureCounts")

# Viewing Files in the Directory
list.files()

# Reading the Count Matrix
countData <- read.delim("./matriz_arabidopsis_2023.txt", header = TRUE, row.names = 1)
head(countData)

#d. Description of Samples
condition <- factor(c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment"))
colData <- data.frame(row.names = colnames(countData), condition)
head(colData)

# Creating a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
dds

#e. Generating a PCA Plot
rld <- rlog(dds, blind = F)
plotPCA(rld, intgroup = "condition") + geom_text(aes(label=name),
                                                 vjust=0.2)


#f. Filtering Genes with Very Low Expression (less than 10 reads)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds

# Performing Differential Expression Analysis
dds <- DESeq(dds)
dds

# Extracting Results
res <- results(dds)
res

# Summary of Differential Gene Expression
summary(res)

# Sorting the Summary List by Adjusted p-adj
res <- res[order(res$padj), ]
head(res)


#g. Extracting Contrasts between Conditions
Treatment_vs_Control <- results(dds, contrast = c("condition", "Treatment", "Control"))
summary(Treatment_vs_Control)

# Obtaining a List of Differentially Expressed Genes with Stricter Filtering
deg <- subset(Treatment_vs_Control, padj < 0.05 & abs(log2FoldChange) > 1)
print(deg)
summary(deg)


#h. Exporting the DEG Table
write.csv(deg, file = "DEG_Treatment_vs_Control_strict.csv")

# Printing and Exporting up-DEGs
up <- subset(deg, log2FoldChange > 1)
print(up)
summary(up)
write.csv(up, file = "Up_Treatment_vs_Control_strict.csv")

# Printing and Exporting down-DEGs
down <- subset(deg, log2FoldChange < (-1))
print(down)
summary(down)
write.csv(down, file = "Down_Treatment_vs_Control_strict.csv")

# Generating a Volcano Plot with EnhancedVolcano
EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'pvalue')

# Setting the Contrast of Interest
plotMA(Treatment_vs_Control, alpha = 0.05, main = "Inoculated with Trichoderma vs Control", xlab = "mean of normalized counts")


# Selecting the Top 20 Differentially Expressed Genes in the Treatment_vs_Control Comparison
res_ordered <- Treatment_vs_Control[order(Treatment_vs_Control$padj), ]
top_genes <- row.names(res_ordered)[1:20]

# Extracting and Normalizing Counts
counts <- counts(dds, normalized = TRUE)
counts_top <- counts[top_genes, ]

# Applying a Logarithmic Transformation to Counts
log_counts_top <- log2(counts_top + 1)

# Creating an Annotation Data Frame Based on Condition Information (colData)
df <- colData

# Displaying the Annotation Data Frame
df

# Generating a Heatmap using the pheatmap library

heatmap_20 <- pheatmap(log_counts_top, annotation = df)











## Protocol 2

### De novo assembly of non-model plant

In this protocol, we will be working with a dataset from Agave sisalana, obtained from a study by Sarwar et al., 2019, where they investigated drought tolerance in this plant. Using these data, we will perform de novo assembly and differential expression analysis. Trinity will be used as the assembler, and DESeq2 will be employed for the differential expression analysis (Haas et al., 2013; Love et al., 2014). To execute this protocol on our personal computer, we have trimmed the original data to only 10,000 reads.

The dataset cited in Pola-Sánchez et al., 2024, for performing the De Novo assembly can be downloaded from this repository:  

  

Library 1:  
   FastqR1: https://github.com/jmvillalobos/RNA-seq-protocol/tree/main#:~:text=Tra_SRR5137658_1_P.fastq.gz  
   FastqR2: https://github.com/jmvillalobos/RNA-seq-protocol/tree/main#:~:text=Tra_SRR5137658_2_P.fastq.gz  

Library 2:
   FastqR1: https://github.com/jmvillalobos/RNA-seq-protocol/tree/main#:~:text=Tra_SRR5137660_1_P.fastq.gz  
   FastqR2: https://github.com/jmvillalobos/RNA-seq-protocol/tree/main#:~:text=Tra_SRR5137660_2_P.fastq.gz  

Library 3:    
   FastqR1: https://github.com/jmvillalobos/RNA-seq-protocol/tree/main#:~:text=Tra_SRR5137663_1_P.fastq.gz  
   FastqR2: https://github.com/jmvillalobos/RNA-seq-protocol/tree/main#:~:text=Tra_SRR5137663_2_P.fastq.gz  

Library 4:   
   FastqR1: https://github.com/jmvillalobos/RNA-seq-protocol/tree/main#:~:text=Con_SRR5137659_1_P.fastq.gz  
   FastqR2: https://github.com/jmvillalobos/RNA-seq-protocol/tree/main#:~:text=Con_SRR5137659_2_P.fastq.gz  

Library 5:   
   FastqR1: https://github.com/jmvillalobos/RNA-seq-protocol/tree/main#:~:text=Con_SRR5137661_1_P.fastq.gz  
   FastqR2: https://github.com/jmvillalobos/RNA-seq-protocol/tree/main#:~:text=Con_SRR5137661_2_P.fastq.gz  

Library 6:  
   FastqR1: https://github.com/jmvillalobos/RNA-seq-protocol/tree/main#:~:text=Con_SRR5137662_1_P.fastq.gz  
   FastqR2: https://github.com/jmvillalobos/RNA-seq-protocol/tree/main#:~:text=Con_SRR5137662_2_P.fastq.gz  
   
