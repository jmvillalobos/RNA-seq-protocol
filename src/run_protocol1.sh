#!/bin/bash

RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

if [ ! -f "/RNA_protocol/quantification_featureCounts/matriz_arabidopsis_2023.txt" ]; then
    wget -P /RNA_protocol/quantification_featureCounts/ https://zenodo.org/records/10576137/files/matriz_arabidopsis_2023.txt

echo -e "${RED}INITIATING PROTOCOL 1${NC}"

echo -e "${BLUE}Running Differential Expression Analysis...${NC}"
(cd /RNA_protocol/quantification_featureCounts/ && Rscript /src/scripts/differential_expression_analysis.R)

echo -e "${RED}PROTOCOL 1 COMPLETED${NC}"
