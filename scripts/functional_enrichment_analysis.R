library("biomartr")
library("clusterProfiler")
library("tidyverse")
library("enrichplot")
suppressPackageStartupMessages(library("org.At.tair.db"))

pdf(file="PlotsFuntionalEnrichment.pdf")

diff_genes <- read_delim(file = "Up_Treatment_vs_Control_strict.csv", delim = ",")

colnames(diff_genes)[1] <- "genes"

diff_genes <- diff_genes[, c("genes", "log2FoldChange")]

write.table(diff_genes, file = "diff_genes_up.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


attributes_to_retrieve <- c("tair_symbol", "entrezgene_id")


result_BM <- biomartr::biomart(
    genes = diff_genes$genes,
    mart = "plants_mart",
    dataset = "athaliana_eg_gene",
    attributes = attributes_to_retrieve,
    filters = "ensembl_gene_id"
)

all_arabidopsis_genes <- read.delim("matriz_arabidopsis_2023.txt", header = TRUE, stringsAsFactors = FALSE)[,1]

all_arabidopsis_genes_annotated <- biomartr::biomart(
    genes = all_arabidopsis_genes,
    mart = "plants_mart",
    dataset = "athaliana_eg_gene",
    attributes = attributes_to_retrieve,
    filters = "ensembl_gene_id"
)

all_arabidopsis_genes_annotated$entrezgene_id <- as.character(all_arabidopsis_genes_annotated$entrezgene_id)

diff_arabidopsis_genes_annotated <- biomartr::biomart(
    genes = diff_genes$genes,
    mart = "plants_mart",
    dataset = "athaliana_eg_gene",
    attributes = attributes_to_retrieve,
    filters = "ensembl_gene_id"
)

ora_analysis_bp <- enrichGO(
    gene = diff_arabidopsis_genes_annotated$entrezgene_id,
    universe = all_arabidopsis_genes_annotated$entrezgene_id,
    OrgDb = org.At.tair.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE,
    pool = FALSE
)

ora_analysis_bp_simplified <- clusterProfiler::simplify(ora_analysis_bp)

write_delim(
    x = as.data.frame(ora_analysis_bp_simplified@result),
    path = "go_results_up.csv",
    delim = ","
)

dotplot(ora_analysis_bp_simplified, showCategory = 10)

barplot(ora_analysis_bp_simplified, showCategory = 10)

ora_analysis_bp <- pairwise_termsim(ora_analysis_bp, method = "JC")

emapplot(ora_analysis_bp, color = "qvalue", showCategory = 15)

