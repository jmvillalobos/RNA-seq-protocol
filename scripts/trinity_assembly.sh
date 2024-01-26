#!/bin/bash

Trinity --seqType fq --samples_file /RNA_protocol/novo_assembly/trinity_analysis/samples.txt --CPU 4 --max_memory 12G --output /RNA_protocol/novo_assembly/trinity_analysis/

mv ../RNA_protocol/novo_assembly/trinity_analysis/trinity_out_dir.Trinity.fasta ../RNA_protocol/novo_assembly/trinity_analysis/Trinity.fasta 
mv ../RNA_protocol/novo_assembly/trinity_analysis/trinity_out_dir.Trinity.fasta.gene_trans_map ../RNA_protocol/novo_assembly/trinity_analysis/Trinity.fasta.gene_trans_map
