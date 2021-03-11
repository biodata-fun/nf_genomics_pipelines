# This image installs BEDTools, Samtools, tabix, bgzip, python, wget, git, etc... 
#parent image
FROM ubuntu:latest

LABEL maintainer="ernestolowy@gmail.com"
LABEL description="Dockerfile used to build an image containing the SAMTools (http://www.htslib.org/) and derived (Tabix, BGZIP) dependencies" 

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get -qq install \
                    wget \
    				build-essential \
					tabix \
				    autoconf \
				    zlib1g-dev \
				    libbz2-dev \
				    liblzma-dev \
				    libhts-dev  \
				    libcurl4-openssl-dev \
				    libssl-dev \
				    && apt-get clean

WORKDIR tmp/

#install Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 && tar -xvf samtools-1.11.tar.bz2 
WORKDIR samtools-1.11
RUN ./configure --without-curses && make && make install
#cleaning
WORKDIR /tmp/
RUN rm -r samtools-1.11/ && rm samtools-1.11.tar.bz2
