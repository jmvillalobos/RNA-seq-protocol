#!/bin/bash

RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

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

echo -e "${RED}INITIATING PROTOCOL 3${NC}"

echo -e "${BLUE}Undertaking de novo assebly with Trinity...${NC}"
bash /src/scripts/trinity_assembly.sh 2>&1 | tee /RNA_protocol/script_logs/trinity_assembly.log

echo -e "${BLUE}Using CD-HIT for redundancy removal of the assembly...${NC}"
bash /src/scripts/cdhit_redundancyRemoval.sh 2>&1 | tee /RNA_protocol/script_logs/cdhit_redundancyRemoval.log

echo -e "${BLUE}Assesing assembly quality...${NC}"
bash /src/scripts/trinityStats_analysis.sh 2>&1 | tee /RNA_protocol/script_logs/trinityStats_analysis.log

echo -e "${BLUE}Performing transcript expression analysis with Salmon...${NC}"
bash /src/scripts/salmon_transcriptExpression.sh 2>&1 | tee /RNA_protocol/script_logs/salmon_transcriptExpression.log

echo -e "${RED}PROTOCOL 3 COMPLETED${NC}"