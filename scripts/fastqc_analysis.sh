#!/bin/bash

for sample in /RNA_protocol/raw_data/*.gz
do
    filename=$(basename ${sample})
    fastqc "/RNA_protocol/raw_data/${filename}" --outdir "/RNA_protocol/quality_raw/"
done