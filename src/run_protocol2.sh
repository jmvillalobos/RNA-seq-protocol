#!/bin/bash

RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${RED}INITIATING PROTOCOL 2${NC}"

echo -e "${BLUE}Performing feature quantification with featureCounts...${NC}"
bash /usr/scripts/featureCounts_quantification.sh 2>&1 | tee /RNA_protocol/script_logs/featureCounts_quantification.log

echo -e "${BLUE}Running Differential Expression Analysis...${NC}"
(cd /RNA_protocol/quantification_featureCounts/ && Rscript /usr/scripts/differential_expression_analysis.R)

echo -e "${BLUE}Running Functional Enrichment Analysis...${NC}"
(cd /RNA_protocol/quantification_featureCounts/ && Rscript /usr/scripts/functional_enrichment_analysis.R)

echo -e "${RED}PROTOCOL 2 COMPLETED${NC}"