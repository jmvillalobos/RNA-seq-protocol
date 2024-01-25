#!/bin/bash

hisat2-build -p 4 /RNA_protocol/genome_arabidopsis/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa /RNA_protocol/index_hisat2/genome

for sample in 204 206 210 212 216 218; 
do
    hisat2 -p 4 -x /RNA_protocol/index_hisat2/genome \
        -1 /RNA_protocol/trimming_data/SRR10207${sample}_P_1.fastq.gz \
        -2 /RNA_protocol/trimming_data/SRR10207${sample}_P_2.fastq.gz \
        -S /RNA_protocol/alignment_hisat2/SRR10207${sample}.sam

    samtools sort /RNA_protocol/alignment_hisat2/SRR10207${sample}.sam -o  /RNA_protocol/alignment_hisat2/SRR10207${sample}_sorted.bam
done

