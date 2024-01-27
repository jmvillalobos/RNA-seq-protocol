#!/bin/bash

RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${RED}INITIATING PROTOCOL 1${NC}"

echo -e "${BLUE}Downloading data...${NC}"
bash /src/scripts/download_data.sh 
# bash /src/scripts/download_data.sh 2>&1 | tee /RNA_protocol/script_logs/download_data.log

echo -e "${BLUE}Performing data quality analysis with FastQC...${NC}"
bash /src/scripts/fastqc_analysis.sh 2>&1 | tee /RNA_protocol/script_logs/fastqc_analysis.log

echo -e "${BLUE}Trimming low quality sequences using Trimmomatic...${NC}"
bash /src/scripts/trimmomatic_filtering.sh 2>&1 | tee /RNA_protocol/script_logs/trimmomatic_filtering.log

echo -e "${BLUE}Aligning reads to reference genome...${NC}"
bash /src/scripts/hisat2_alignment.sh 2>&1 | tee /RNA_protocol/script_logs/hisat2_alignment.log

echo -e "${RED}PROTOCOL 1 COMPLETED${NC}"
