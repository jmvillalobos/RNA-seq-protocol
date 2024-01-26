#!/bin/bash

data_path="/RNA_protocol/alignment_hisat2/"

for file in ${data_path}/*.bam; do
    name=$(basename ${file} .bam)
    featureCounts -p --countReadPairs -t exon -g gene_id \
        -a /RNA_protocol/genome_arabidopsis/Arabidopsis_thaliana.TAIR10.57.gtf \
        -o /RNA_protocol/quantification_featureCounts/featureCounts_output/${name}.txt /RNA_protocol/alignment_hisat2/${name}.bam
done

ls -1 /RNA_protocol/quantification_featureCounts/featureCounts_output/*.txt | parallel 'cat {} | sed '1d' | cut -f7 {} > /RNA_protocol/quantification_featureCounts/{/.}_clean.txt'
ls -1 /RNA_protocol/quantification_featureCounts/featureCounts_output/*.txt | head -1 | xargs cut -f1 > /RNA_protocol/quantification_featureCounts/genes.txt
paste /RNA_protocol/quantification_featureCounts/genes.txt /RNA_protocol/quantification_featureCounts/*clean.txt > /RNA_protocol/quantification_featureCounts/matriz_arabidopsis_2023.txt