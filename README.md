# RNA-seq-protocol
 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10537097.svg)](https://zenodo.org/DOI/10.5281/zenodo.10537096)


 This is the protocol for conducting RNAseq analysis in plants. This repository has been created to make it easier to follow the protocol presented in the Current Protocols journal.

  To streamline the execution of the scripts within this repository, we recommend using Docker. Docker provides a consistent environment for running applications, ensuring that your scripts run smoothly across different platforms. If you haven't installed Docker yet, you can download it based on your operating system:

- **Linux Users:** [Install Docker for Linux](https://docs.docker.com/desktop/install/linux-install/)
- **Windows Users:** [Install Docker for Windows](https://docs.docker.com/desktop/install/windows-install/)
- **MacOS Users:** [Install Docker for Mac](https://docs.docker.com/desktop/install/mac-install/)

Once Docker is installed, you'll be able to effortlessly run the scripts in this repository and benefit from a reproducible and isolated environment. Follow the steps in our [Running](#run) section to reproduce our results.


## Running {#run}


<!-- ![Docker](https://img.shields.io/badge/docker-%230db7ed.svg?style=for-the-badge&logo=docker&logoColor=white) -->



## Protocol 1


**Step 1: Data Download**
Download raw RNA-seq data from the European Nucleotide Archive (ENA) related to the study of Trichoderma atroviride's Nox genes during interaction with Arabidopsis thaliana.

**Step 3: Quality Check on Raw Reads with FastQC**

Evaluate the quality of raw data using FastQC, analyzing metrics such as per base sequence quality.

**Step 4: Read Trimming with Trimmomatic**

Trim raw reads using Trimmomatic to remove low-quality bases, adapter sequences, and short reads.

**Step 5: Alignment of Reads**

- Create an index for the reference genome.
- Align trimmed reads to the Arabidopsis thaliana reference genome using HISAT2.

**Step 6: Quantification of Mapped Reads**

- Download the GTF annotation file for Arabidopsis thaliana.
- Quantify the number of mapped reads using featureCounts.
- Create a count matrix for downstream analysis.

**Step 7: Differential Expression Analysis**

- Perfor differential expression analysis using DESeq2 in R.
- Filter lowly expressed genes.
- Generate PCA plots for quality assessment.
- Exporte results, create volcano plots, and identify differentially expressed genes (DEGs).

## Protocol 2

**Functional Enrichment Analysis**




<p align="center">
<img src="imgs/protocol1and2.png" alt="GitHub Logo" width="400" height="490">
</p>

## Protocol 3
**De novo Assembly with Trinity**
- Executed de novo transcriptome assembly using Trinity.
- Resulted in `Trinity.fasta` representing assembled transcripts.

**Step 2: Evaluating the Assembly**

- Count Assembled Transcripts

- Used a command to count the number of assembled transcripts and obtain Comprehensive Assembly Statistics

- Included N50, GC content, and contig lengths.

**Step 3: Redundancy Removal**

- Applied CD-HIT to remove redundancy from the transcriptome assembly.
- Reduced the number of contigs.

**Step 4: Transcript Expression Quantification with Salmon**

- Used Salmon to estimate expression values of transcripts.
- Created directories for each replicate and inspected the quantification results.


<p align="center">
<img src="imgs/protocol3.png" alt="GitHub Logo" width="400" height="420">
</p>
### Raw Data for the analysis:

In the 7th step of the protocol outlined in Current Protocols, we provide detailed instructions for conducting a differential expression analysis. To reach this stage, it is imperative to download the Fastq files and complete the entire analysis, culminating in the contingency table containing counts for each gene. However, to simplify the initiation of hands-on experience with statistical analyses in R, we have deposited the final count matrix here for your convenience.

Matrix: [Matrix_arabidopsis_2023](https://raw.githubusercontent.com/jmvillalobos/RNA-seq-protocol/main/matriz_arabidopsis_2023.txt)


#Supplementary 2

#Differential gene expression analysis 

#a.Library Installation
#Installing DESeq2
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

#Installing ggplot2
install.packages("ggplot2")

#Installing EnhancedVolcano
BiocManager::install("EnhancedVolcano")

#Installing pheatmap
install.packages("pheatmap")




#b. Loading the library for DEG
library(DESeq2)
#Loading the library for plots
library("ggplot2")
#Loading the library for volcano plot
library("EnhancedVolcano")
#Loading the library for heatmap
library("pheatmap")




#c.Setting the Path to the featureCounts Count Matrix
setwd("C:/project_2023/quantification_featureCounts")

#Viewing Files in the Directory
list.files()

#Reading the Count Matrix
countData <- read.delim("./matriz_arabidopsis_2023.txt", header = TRUE, row.names = 1)
head(countData)




#d. Description of Samples
condition <- factor(c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment"))
colData <- data.frame(row.names = colnames(countData), condition)
head(colData)

#Creating a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
dds




#e. Generating a PCA Plot
rld <- rlog(dds, blind = F)
plotPCA(rld, intgroup = "condition") + geom_text(aes(label=name),

                                                 vjust=0.2)


#f. Filtering Genes with Very Low Expression (less than 10 reads)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds

#Performing Differential Expression Analysis
dds <- DESeq(dds)
dds

#Extracting Results
res <- results(dds)
res

#Summary of Differential Gene Expression
summary(res)

#Sorting the Summary List by Adjusted p-adj
res <- res[order(res$padj), ]
head(res)




#g. Extracting Contrasts between Conditions
Treatment_vs_Control <- results(dds, contrast = c("condition", "Treatment", "Control"))
summary(Treatment_vs_Control)

#Obtaining a List of Differentially Expressed Genes with Stricter Filtering
deg <- subset(Treatment_vs_Control, padj < 0.05 & abs(log2FoldChange) > 1)
print(deg)
summary(deg)




#h. Exporting the DEG Table
write.csv(deg, file = "DEG_Treatment_vs_Control_strict.csv")

#Printing and Exporting up-DEGs
up <- subset(deg, log2FoldChange > 1)
print(up)
summary(up)
write.csv(up, file = "Up_Treatment_vs_Control_strict.csv")

#Printing and Exporting down-DEGs
down <- subset(deg, log2FoldChange < (-1))
print(down)
summary(down)
write.csv(down, file = "Down_Treatment_vs_Control_strict.csv")

#Generating a Volcano Plot with EnhancedVolcano
EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'pvalue')

#Setting the Contrast of Interest
plotMA(Treatment_vs_Control, alpha = 0.05, main = "Inoculated with Trichoderma vs Control", xlab = "mean of normalized counts")


#Selecting the Top 20 Differentially Expressed Genes in the Treatment_vs_Control Comparison
res_ordered <- Treatment_vs_Control[order(Treatment_vs_Control$padj), ]
top_genes <- row.names(res_ordered)[1:20]

#Extracting and Normalizing Counts
counts <- counts(dds, normalized = TRUE)
counts_top <- counts[top_genes, ]

#Applying a Logarithmic Transformation to Counts
log_counts_top <- log2(counts_top + 1)

#Creating an Annotation Data Frame Based on Condition Information (colData)
df <- colData

#Displaying the Annotation Data Frame
df

#Generating a Heatmap using the pheatmap library

heatmap_20 <- pheatmap(log_counts_top, annotation = df)










## Protocol 2

### Enrichment analysis 

#Supplementary 3
#Gene ontology enrichment analysis

#a.Installation of the package to use in the R session.
#Installing clusterProfiler
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

#Installing biomartr 1.0.7 from CRAN
install.packages("biomartr", dependencies = TRUE)


#Installing tidyverse
install.packages("tidyverse")

#Installing enrichplot 
BiocManager::install("enrichplot")




#b.Load required libraries
library("biomartr")
library("clusterProfiler")
library("tidyverse")
library("enrichplot")
suppressPackageStartupMessages(library("org.At.tair.db"))




#c. Load the DESeq2 table into the R session.

#c.1. Setting the path to the DESeq2 table.
setwd("C:/project_2023/quantification_featureCounts")
list.files()

#c.2. Read the DESeq2 table (up or down).
diff_genes <- read_delim(file = "Up_Treatment_vs_Control_strict.csv", delim = ",")

#c.3. Assign names to the first column.
colnames(diff_genes)[1] <- "genes"

#c.4. Create a new table with the columns of interest.
diff_genes <- diff_genes[, c("genes", "log2FoldChange")]

#c.5. Save the new table to a file.
write.table(diff_genes, file = "diff_genes_up.tsv", sep = "\t", row.names = FALSE, quote = FALSE)




#c. DEG annotation with Ensembl and biomartr.

#c.1. We retrieved the athaliana_eg_gene dataset from the TAIR10 genome annotation version of the plant_mart mart.
biomartr::organismBM(organism = "Arabidopsis thaliana")

#c.2. We recovered the correspondence between our Arabidopsis differentially expressed gene identifiers (Ensembl gene id) with tair_symbol and entrezgene_id.

attributes_to_retrieve <- c("tair_symbol", "entrezgene_id")

result_BM <- biomartr::biomart(
  genes = diff_genes$genes,
  mart = "plants_mart",
  dataset = "athaliana_eg_gene",
  attributes = attributes_to_retrieve,
  filters = "ensembl_gene_id"
)

head(result_BM)




#d. Over-Representation Analysis (ORA) with ClusterProfiler.

#d.1. Change working directory.
setwd("C:/project_2023/quantification_featureCounts")

#d.2. Read all Arabidopsis genes from a file (universe).
all_arabidopsis_genes <- read.delim("matriz_arabidopsis_2023.txt", header = TRUE, stringsAsFactors = FALSE)[,1]

#d.3. We annotated our gene universe by querying the Ensembl API using the biomartr::biomart function with the 'athaliana_eg_gene' dataset. The resulting table includes the original ensembl_gene_id, tair_symbol, and entrezgene_id for each gene in our dataset.
all_arabidopsis_genes_annotated <- biomartr::biomart(
  genes = all_arabidopsis_genes,
  mart = "plants_mart",
  dataset = "athaliana_eg_gene",
  attributes = attributes_to_retrieve,
  filters = "ensembl_gene_id"
)

head(all_arabidopsis_genes_annotated)


#d.4. For compatibility with the enrichGO function, all_arabidopsis_genes_annotated genes must be characters, not integers, so their conversion is necessary.
all_arabidopsis_genes_annotated$entrezgene_id <- as.character(
  all_arabidopsis_genes_annotated$entrezgene_id
)


#d.5. Perform ORA for Gene Ontology Biological Process class.
ora_analysis_bp <- enrichGO(
  gene = result_BM$entrezgene_id,
  universe = all_arabidopsis_genes_annotated$entrezgene_id,
  OrgDb = org.At.tair.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE,
  pool = FALSE
)



#d.6. Simplify the ORA results.
ora_analysis_bp_simplified <- clusterProfiler::simplify(ora_analysis_bp)




#e. Exporting and plotting figures.

#e.1. Write simplified results to a CSV file.
write_delim(
  x = as.data.frame(ora_analysis_bp_simplified@result),
  path = "go_results_up.csv",
  delim = ","
)

#e.2. Plot 1 - Dotplot.
dotplot(ora_analysis_bp_simplified, showCategory = 10)


#e.3. Plot 2 - Barplot.
barplot(ora_analysis_bp_simplified, showCategory = 10)

#e.4 Plot 3 - Emap Plot.
ora_analysis_bp <- pairwise_termsim(ora_analysis_bp, method = "JC")
emapplot(ora_analysis_bp, color = "qvalue", showCategory = 15)










## Protocol 3

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
   
