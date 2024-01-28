#!/bin/bash

wget /RNA_protocol/raw_data/ https://zenodo.org/records/10576137/files/raw_data1.zip
wget /RNA_protocol/raw_data/ https://zenodo.org/records/10576137/files/raw_data2.zip
wget /RNA_protocol/raw_data/ https://zenodo.org/records/10576137/files/raw_data3.zip

unzip /RNA_protocol/raw_data/raw_data1.zip -d /RNA_protocol/raw_data/
unzip /RNA_protocol/raw_data/raw_data2.zip -d /RNA_protocol/raw_data/
unzip /RNA_protocol/raw_data/raw_data3.zip -d /RNA_protocol/raw_data/
rm /RNA_protocol/raw_data/raw_data1.zip
rm /RNA_protocol/raw_data/raw_data2.zip
rm /RNA_protocol/raw_data/raw_data3.zip

wget -P /RNA_protocol/genome_arabidopsis/ -nc https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-57/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gunzip /RNA_protocol/genome_arabidopsis/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

wget -P /RNA_protocol/genome_arabidopsis/ -nc https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-57/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.57.gtf.gz
gunzip /RNA_protocol/genome_arabidopsis/Arabidopsis_thaliana.TAIR10.57.gtf.gz

wget -P /RNA_protocol/novo_assembly/trimming_data_agave/ -nc https://zenodo.org/records/10576137/files/Con_SRR5137659_1_P.fastq.gz
wget -P /RNA_protocol/novo_assembly/trimming_data_agave/ -nc https://zenodo.org/records/10576137/files/Con_SRR5137659_2_P.fastq.gz
wget -P /RNA_protocol/novo_assembly/trimming_data_agave/ -nc https://zenodo.org/records/10576137/files/Con_SRR5137661_1_P.fastq.gz
wget -P /RNA_protocol/novo_assembly/trimming_data_agave/ -nc https://zenodo.org/records/10576137/files/Con_SRR5137661_2_P.fastq.gz
wget -P /RNA_protocol/novo_assembly/trimming_data_agave/ -nc https://zenodo.org/records/10576137/files/Con_SRR5137662_1_P.fastq.gz
wget -P /RNA_protocol/novo_assembly/trimming_data_agave/ -nc https://zenodo.org/records/10576137/files/Con_SRR5137662_2_P.fastq.gz
wget -P /RNA_protocol/novo_assembly/trimming_data_agave/ -nc https://zenodo.org/records/10576137/files/Tra_SRR5137658_1_P.fastq.gz
wget -P /RNA_protocol/novo_assembly/trimming_data_agave/ -nc https://zenodo.org/records/10576137/files/Tra_SRR5137658_2_P.fastq.gz
wget -P /RNA_protocol/novo_assembly/trimming_data_agave/ -nc https://zenodo.org/records/10576137/files/Tra_SRR5137660_1_P.fastq.gz
wget -P /RNA_protocol/novo_assembly/trimming_data_agave/ -nc https://zenodo.org/records/10576137/files/Tra_SRR5137660_2_P.fastq.gz
wget -P /RNA_protocol/novo_assembly/trimming_data_agave/ -nc https://zenodo.org/records/10576137/files/Tra_SRR5137663_1_P.fastq.gz
wget -P /RNA_protocol/novo_assembly/trimming_data_agave/ -nc https://zenodo.org/records/10576137/files/Tra_SRR5137663_2_P.fastq.gz