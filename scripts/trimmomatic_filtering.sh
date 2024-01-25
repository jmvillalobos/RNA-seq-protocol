for sample in 204 206 210 212 216 218; 
do
    trimmomatic PE -phred33 \
        /RNA_protocol/raw_data/SRR10207${sample}_1.fastq.gz \
        /RNA_protocol/raw_data/SRR10207${sample}_2.fastq.gz \
        /RNA_protocol/trimming_data/SRR10207${sample}_P_1.fastq.gz /RNA_protocol/trimming_data/SRR10207${sample}_U_1.fastq.gz \
        /RNA_protocol/trimming_data/SRR10207${sample}_P_2.fastq.gz /RNA_protocol/trimming_data/SRR10207${sample}_U_2.fastq.gz \
        ILLUMINACLIP:/usr/local/bin/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
        SLIDINGWINDOW:4:15 MINLEN:36
done