# This image installs the basic Ubuntu image without more dependencies 
FROM ubuntu:latest

LABEL maintainer="ernestolowy@gmail.com"
LABEL description="Dockerfile used to build an Ubuntu-based image used by different workflows" 

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get -qq install \
    				build-essential \
				    autoconf \
				    && apt-get clean
