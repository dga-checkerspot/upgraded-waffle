FROM ubuntu:18.04
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN apt-get update

RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*


RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 
RUN conda --version

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

RUN conda install -c conda-forge -y awscli=1.22.77
RUN conda install -c bioconda sra-tools=2.11.0
RUN conda install -c bioconda bcftools=1.15.1
RUN conda install -c bioconda bwa=0.7.17
RUN conda install -c bioconda samtools=1.15.1
RUN conda install -c bioconda trinity
RUN conda install -c bioconda cutadapt
RUN conda install -c bioconda bbmap
RUN conda install -c bioconda fastq-pair
RUN conda install -c bioconda spades
RUN conda install -c bioconda seqtk
