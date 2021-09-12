FROM alpine:3.14

LABEL maintainer="ernestolowy@gmail.com"
LABEL description="Dockerfile used to build BCFTools" 

RUN apk add --no-cache git build-base zlib-dev bzip2-dev xz-dev curl-dev

WORKDIR tmp/

# install BCFTools
RUN git clone --recurse-submodules git://github.com/samtools/htslib.git && git clone git://github.com/samtools/bcftools.git
WORKDIR bcftools
RUN make && make install
WORKDIR /tmp/
RUN rm -rf htslib && rm -rf bcftools