# RNA-seq-protocol
 This is the protocol for conducting RNAseq analysis in plants. This repository has been created to make it easier to follow the protocol presented in the Current Protocols journal.
## Protocol 1

### Raw Data for the analysis:

In the 7th step of the protocol outlined in Current Protocols, we provide detailed instructions for conducting a differential expression analysis. To reach this stage, it is imperative to download the Fastq files and complete the entire analysis, culminating in the contingency table containing counts for each gene. However, to simplify the initiation of hands-on experience with statistical analyses in R, we have deposited the final count matrix here for your convenience.

Matrix: [Matrix_arabidopsis_2023](https://raw.githubusercontent.com/jmvillalobos/RNA-seq-protocol/main/matriz_arabidopsis_2023.txt)


## Protocol 2

### De novo assembly of non-model plant

In this protocol, we will be working with a dataset from Agave sisalana, obtained from a study by Sarwar et al., 2019, where they investigated drought tolerance in this plant. Using these data, we will perform de novo assembly and differential expression analysis. Trinity will be used as the assembler, and DESeq2 will be employed for the differential expression analysis (Haas et al., 2013; Love et al., 2014). To execute this protocol on our personal computer, we have trimmed the original data to only 10,000 reads.

The dataset cited in Pola-Sánchez et al., 2024, for performing the De Novo assembly can be downloaded from this repository:

Sample 1:
   FastqR1: https://github.com/jmvillalobos/RNA-seq-protocol/tree/main#:~:text=Tra_SRR5137658_1_P.fastq.gz
   FastqR2: https://github.com/jmvillalobos/RNA-seq-protocol/tree/main#:~:text=Tra_SRR5137658_2_P.fastq.gz
