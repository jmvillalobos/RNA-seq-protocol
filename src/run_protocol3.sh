#!/bin/bash

RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

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