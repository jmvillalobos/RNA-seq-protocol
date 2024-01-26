# Use the official Ubuntu base image
FROM ubuntu:22.04 as intermediate

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y zip wget

# TRIMMOMATIC v0.39
RUN apt-get install -y openjdk-8-jre-headless
RUN wget https://github.com/usadellab/Trimmomatic/files/5854859/Trimmomatic-0.39.zip
RUN unzip Trimmomatic-0.39.zip
RUN mv Trimmomatic-0.39 /usr/local/bin/
RUN touch /usr/local/bin/trimmomatic
RUN echo '#!/bin/bash' > /usr/local/bin/trimmomatic
RUN echo 'java -jar /usr/local/bin/Trimmomatic-0.39/trimmomatic-0.39.jar "$@"' >> /usr/local/bin/trimmomatic
RUN chmod +x /usr/local/bin/trimmomatic
RUN rm -r Trimmomatic-0.39.zip


# STAR v2.7.10a
RUN wget --no-check-certificate https://github.com/alexdobin/STAR/archive/2.7.10a.zip
RUN unzip 2.7.10a.zip
RUN apt-get install make
RUN apt-get update
RUN apt-get -y install build-essential
RUN apt-get -y install zlib1g-dev
# RUN sed -i 's/arrIn\[ibit\]/arrIn.at(ibit);/g' STAR-2.7.10a/source/SoloCommon.h
RUN sed -i '1i #include <array>' STAR-2.7.10a/source/SoloCommon.h
RUN (cd STAR-2.7.10a/source/ && make)
RUN mv STAR-2.7.10a/ /usr/local/bin/
RUN touch /usr/local/bin/STAR
RUN echo '#!/bin/bash' > /usr/local/bin/STAR
RUN echo '/usr/local/bin/STAR-2.7.10a/source/./STAR' >> /usr/local/bin/STAR
RUN chmod +x /usr/local/bin/STAR
RUN rm 2.7.10a.zip

# FastQC v11.9
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
RUN unzip fastqc_v0.11.9.zip
RUN mv FastQC/ /usr/local/bin/
RUN touch /usr/local/bin/fastqc
RUN echo '#!/bin/bash' > /usr/local/bin/fastqc
RUN echo 'perl /usr/local/bin/FastQC/fastqc "$@"' >> /usr/local/bin/fastqc
RUN chmod +x /usr/local/bin/fastqc
RUN rm fastqc_v0.11.9.zip

#Kallisto v0.46.1
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz
RUN tar -xvf kallisto_linux-v0.46.1.tar.gz
RUN mv kallisto/ /usr/local/bin/kallisto_folder
RUN touch /usr/local/bin/kallisto
RUN echo '#!/bin/bash' > /usr/local/bin/kallisto
RUN echo '/usr/local/bin/kallisto_folder/./kallisto "$@"' >> /usr/local/bin/kallisto
RUN chmod +x /usr/local/bin/kallisto
RUN rm kallisto_linux-v0.46.1.tar.gz

#MultiQC v1.12
RUN apt-get install -y python3-pip
RUN pip3 install --upgrade pip
RUN pip3 install multiqc==1.12 --break-system-packages

#Salmon v1.4.0
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v1.4.0/salmon-1.4.0_linux_x86_64.tar.gz
RUN tar -xvzf salmon-1.4.0_linux_x86_64.tar.gz
RUN mv salmon-latest_linux_x86_64/ /usr/local/bin/
RUN touch /usr/local/bin/salmon
RUN echo '#!/bin/bash' > /usr/local/bin/salmon
RUN echo '/usr/local/bin/salmon-latest_linux_x86_64/bin/./salmon "$@"' >> /usr/local/bin/salmon
RUN chmod +x /usr/local/bin/salmon
RUN rm salmon-1.4.0_linux_x86_64.tar.gz

#featureCounts/subread v2.0.3
RUN wget https://sourceforge.net/projects/subread/files/subread-2.0.3/subread-2.0.3-Linux-x86_64.tar.gz/download
RUN tar -xvf download
RUN mv subread-2.0.3-Linux-x86_64/ /usr/local/bin/
RUN touch /usr/local/bin/featureCounts
RUN echo '#!/bin/bash' > /usr/local/bin/featureCounts
RUN echo '/usr/local/bin/subread-2.0.3-Linux-x86_64/bin/./featureCounts "$@"' >> /usr/local/bin/featureCounts
RUN chmod +x /usr/local/bin/featureCounts
RUN rm download

#htseq-count v1.99.2
RUN pip3 install cython --break-system-packages
RUN pip3 install HTSeq==1.99.2 --break-system-packages

#Trinity v2.13.2
RUN apt-get -y install cmake
RUN apt-get -y install m4
RUN apt-get -y install libbz2-dev
RUN apt-get -y install liblzma-dev
RUN wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/Trinity-v2.13.2/trinityrnaseq-v2.13.2.FULL.tar.gz
RUN tar -zxvf trinityrnaseq-v2.13.2.FULL.tar.gz
RUN wget https://ftp.gnu.org/gnu/autoconf/autoconf-2.69.tar.gz
RUN tar -zxvf autoconf-2.69.tar.gz
RUN (cd autoconf-2.69 && ./configure && make && make install)
RUN (cd trinityrnaseq-v2.13.2 && make)
RUN mv trinityrnaseq-v2.13.2/ usr/local/bin/
RUN touch /usr/local/bin/Trinity
RUN echo '#!/bin/bash' > /usr/local/bin/Trinity
RUN echo '/usr/local/bin/trinityrnaseq-v2.13.2/./Trinity "$@"' >> /usr/local/bin/Trinity
RUN chmod +x /usr/local/bin/Trinity
RUN rm trinityrnaseq-v2.13.2.FULL.tar.gz
RUN rm autoconf-2.69.tar.gz

#SOAPdenovo-Trans v1.0.5
RUN wget https://github.com/aquaskyline/SOAPdenovo-Trans/archive/refs/tags/1.0.5.tar.gz
RUN tar -zxvf 1.0.5.tar.gz
RUN (cd SOAPdenovo-Trans-1.0.5/src && make)
RUN mv SOAPdenovo-Trans-1.0.5/ usr/local/bin/
RUN touch /usr/local/bin/SOAPdenovo-Trans
RUN echo '#!/bin/bash' > /usr/local/bin/SOAPdenovo-Trans
RUN echo '/usr/local/bin/SOAPdenovo-Trans-1.0.5/./SOAPdenovo-Trans-31mer "$@"' >> /usr/local/bin/SOAPdenovo-Trans
RUN chmod +x /usr/local/bin/SOAPdenovo-Trans
RUN rm 1.0.5.tar.gz


# Samtools 1.9
RUN apt install -y libbz2-dev
RUN apt-get install -y libncurses5-dev libncursesw5-dev liblzma-dev
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
RUN tar -vxjf samtools-1.9.tar.bz2
RUN (cd samtools-1.9 && make)
RUN mv samtools-1.9 usr/local/bin/
RUN touch /usr/local/bin/samtools
RUN echo '#!/bin/bash' > /usr/local/bin/samtools
RUN echo '/usr/local/bin/samtools-1.9/./samtools "$@"' >> /usr/local/bin/samtools
RUN chmod +x /usr/local/bin/samtools
RUN rm samtools-1.9.tar.bz2

#HISAT2 v2.2.1
RUN wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download
RUN unzip download
RUN mv hisat2-2.2.1/ /usr/local/bin/
RUN touch /usr/local/bin/hisat2
RUN echo '#!/bin/bash' > /usr/local/bin/hisat2
RUN echo 'perl /usr/local/bin/hisat2-2.2.1/hisat2 "$@"' >> /usr/local/bin/hisat2
RUN chmod +x /usr/local/bin/hisat2
RUN touch /usr/local/bin/hisat2-build
RUN echo '#!/bin/bash' > /usr/local/bin/hisat2-build
RUN echo 'python3 /usr/local/bin/hisat2-2.2.1/hisat2-build "$@"' >> /usr/local/bin/hisat2-build
RUN chmod +x /usr/local/bin/hisat2-build
RUN rm download


FROM rocker/r-base:4.3.0

RUN apt-get update && \
    apt-get install -y \
        libcurl4-openssl-dev \
        libssl-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libxml2-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev


RUN R -e "install.packages(c('RCurl', 'BiocManager', 'RColorBrewer', 'pheatmap', 'tidyverse', 'ggplot2', 'dplyr'), dependencies=TRUE)" && \
    R -e "BiocManager::install(c('topGO', 'DESeq2', 'edgeR', 'PCAtools', 'EnhancedVolcano', 'clusterProfiler', 'enrichplot', 'biomartr', 'org.At.tair.db'), dependencies=TRUE)" 


RUN apt -y install default-jdk
RUN apt install -y moreutils
RUN apt install -y parallel
RUN apt-get install -y watch
RUN apt-get install -y coreutils
RUN apt install -y cd-hit
RUN apt-get install -y jellyfish
RUN apt-get install -y bowtie2
RUN apt-get install -y python3-pip
RUN pip3 install numpy  --break-system-packages

COPY --from=intermediate /usr/local/bin /usr/local/bin

# Set the working directory
WORKDIR /RNA_protocol


# Set the entry point for the container
ENTRYPOINT ["/bin/bash"]
