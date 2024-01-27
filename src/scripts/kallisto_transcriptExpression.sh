#!/bin/bash

(cd /RNA_protocol/novo_assembly/trinity_analysis/ && perl /usr/local/bin/trinityrnaseq-v2.13.2/util/align_and_estimate_abundance.pl --seqType fq --samples_file /RNA_protocol/novo_assembly/trinity_analysis/samples.txt --transcripts /RNA_protocol/novo_assembly/trinity_analysis/Trinity_90.fasta --est_method kallisto --trinity_mode --prep_reference)

(cd /RNA_protocol/novo_assembly/trinity_analysis/ && find Tra_* Con_* -name "quant.sf" | tee quant_files.list)

(cd /RNA_protocol/novo_assembly/trinity_analysis/ && perl /usr/local/bin/trinityrnaseq-v2.13.2/util/abundance_estimates_to_matrix.pl --est_method salmon --out_prefix Trinity --name_sample_by_basedir --quant_files quant_files.list --gene_trans_map Trinity_90.fasta.gene_trans_map)