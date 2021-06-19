# This image installs BEDTools, Samtools, tabix, bgzip, python, wget, git, etc... 
#parent image
FROM ubuntu:latest

LABEL maintainer="ernestolowy@gmail.com"
LABEL description="Dockerfile used to build an image containing the dependencies used in the normalization.nf nextflow workflow: BCFTools" 

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get -qq install \
                    git \
					wget \
    				build-essential \
				    autoconf \
				    zlib1g-dev \
				    libbz2-dev \
				    liblzma-dev \
				    libhts-dev  \
				    libcurl4-openssl-dev \
				    libssl-dev \
				    && apt-get clean

WORKDIR tmp/

# install BCFTools
RUN git clone --recurse-submodules git://github.com/samtools/htslib.git && git clone git://github.com/samtools/bcftools.git
WORKDIR bcftools
RUN make && make install
WORKDIR tmp/
RUN rm -rf htslib && rm -rf bcftools
