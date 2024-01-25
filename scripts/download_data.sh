#!/bin/bash

# wget -P /RNA_protocol/raw_data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR102/004/SRR10207204/SRR10207204_1.fastq.gz
# wget -P /RNA_protocol/raw_data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR102/004/SRR10207204/SRR10207204_2.fastq.gz
# wget -P /RNA_protocol/raw_data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR102/006/SRR10207206/SRR10207206_1.fastq.gz
# wget -P /RNA_protocol/raw_data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR102/006/SRR10207206/SRR10207206_2.fastq.gz
# wget -P /RNA_protocol/raw_data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR102/010/SRR10207210/SRR10207210_1.fastq.gz
# wget -P /RNA_protocol/raw_data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR102/010/SRR10207210/SRR10207210_2.fastq.gz
# wget -P /RNA_protocol/raw_data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR102/012/SRR10207212/SRR10207212_1.fastq.gz
# wget -P /RNA_protocol/raw_data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR102/012/SRR10207212/SRR10207212_2.fastq.gz
# wget -P /RNA_protocol/raw_data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR102/016/SRR10207216/SRR10207216_1.fastq.gz
# wget -P /RNA_protocol/raw_data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR102/016/SRR10207216/SRR10207216_2.fastq.gz
# wget -P /RNA_protocol/raw_data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR102/018/SRR10207218/SRR10207218_1.fastq.gz
# wget -P /RNA_protocol/raw_data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR102/018/SRR10207218/SRR10207218_2.fastq.gz

# wget -P /RNA_protocol/genome_arabidopsis/ https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-57/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
# gunzip /RNA_protocol/genome_arabidopsis/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

# wget -P /RNA_protocol/genome_arabidopsis/ https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-57/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.57.gtf.gz
# gunzip /RNA_protocol/genome_arabidopsis/Arabidopsis_thaliana.TAIR10.57.gtf.gz

wget -P /RNA_protocol/novo_assembly/trimming_data_agave/ https://zenodo.org/api/records/10537097/files-archive
mv /RNA_protocol/novo_assembly/trimming_data_agave/files-archive /RNA_protocol/novo_assembly/trimming_data_agave/files-archive.zip
unzip /RNA_protocol/novo_assembly/trimming_data_agave/files-archive.zip -d /RNA_protocol/novo_assembly/trimming_data_agave/
rm /RNA_protocol/novo_assembly/trimming_data_agave/files-archive.zip