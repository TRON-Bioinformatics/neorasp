FROM rocker/r-ver:4.4
MAINTAINER johannes.hausmann@tron-mainz.de

ARG TARGETPLATFORM
ENV LANG=EN.UTF-8 LC_ALL=EN.UTF-8

RUN apt-get update && apt-get install -y software-properties-common && \
    apt-get update && apt-get install -y \
        build-essential \
        bzip2 \
        curl \
        git \
        tini \
        libbz2-dev \
        libcurl3-dev \
        liblzma-dev \
        libncurses5-dev \
        libssl-dev \
        unzip \
        wget \
        zlib1g-dev \
        libdeflate-dev \
        libisal-dev \
        libc6-dev \
        python3 \
        python3-pip \
        python-setuptools \
        python3-distutils \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/

RUN ln -sf /usr/bin/python3 /usr/bin/python

ENV SAMTOOLS_VERSION=1.20

RUN SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
   cd $SRC && \
   wget $SAMTOOLS_URL && \
   tar xvf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
   cd samtools-${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION} && ./configure && make && make install && \
   cd ../ && ./configure --without-curses && make && make install && cd $SRC && rm -r samtools-${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}


ENV BOWTIE_VERSION=2.5.4
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v${BOWTIE_VERSION}/bowtie2-${BOWTIE_VERSION}-linux-x86_64.zip -O bowtie2-${BOWTIE_VERSION}-linux-x86_64.zip && \
    unzip bowtie2-${BOWTIE_VERSION}-linux-x86_64.zip && \
    mv bowtie2-${BOWTIE_VERSION}-linux-x86_64/bowtie2* $BIN && \
    rm *.zip && \
    rm -r bowtie2-${BOWTIE_VERSION}-linux-x86_64

ENV STAR_VERSION=2.7.11b
RUN STAR_URL="https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz" &&\
    wget -P $SRC $STAR_URL &&\
    tar -xvf $SRC/${STAR_VERSION}.tar.gz -C $SRC && \
    mv $SRC/STAR-${STAR_VERSION}/bin/Linux_x86_64_static/STAR /usr/local/bin && rm -r $SRC/STAR-${STAR_VERSION}

RUN pip install numpy==1.24.4 build twist pysam==0.21.0 matplotlib==3.5.2 biopython==1.81 bp_quant=0.6.0 snakemake=8.25.2


