#parent image
FROM ubuntu:latest

LABEL maintainer="ernestolowy@gmail.com"
LABEL description="Dockerfile used to build an image used by the Nextflow workflow named cov_per_region.nf Check the nexftlow workflow at https://github.com/biodata-fun/nf_genomics_pipelines/blob/main/workflows/cov_per_region.nf to understand how it works" 

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get -qq install git \
                    wget \
    				build-essential \
				    autoconf \
				    zlib1g-dev \
				    libbz2-dev \
				    liblzma-dev \
				    libhts-dev  \
				    python3 \
				    python3-pip \
				    libcurl4-openssl-dev \
				    libssl-dev \
				    && apt-get clean

WORKDIR tmp/

#prepare Python
RUN ln -s /usr/bin/python3 /usr/bin/python

#install Bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary && mv bedtools.static.binary bedtools && chmod a+x bedtools && mv bedtools /bin/

#install Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 && tar -xvf samtools-1.11.tar.bz2 
WORKDIR samtools-1.11
RUN ./configure --without-curses && make && make install
#cleaning
WORKDIR /tmp/
RUN rm -r samtools-1.11/ && rm samtools-1.11.tar.bz2