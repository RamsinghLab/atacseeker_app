##
## This Docker file builds an image for running the ATACseeker pipeline. 
##

FROM ubuntu:14.04
MAINTAINER Asif Zubair <asif.zubair@gmail.com>
ENV DEBIAN_FRONTEND noninteractive

## Getting ready to install R
RUN echo 'deb http://cran.stat.ucla.edu//bin/linux/ubuntu trusty/' | tee -a /etc/apt/sources.list.d/r.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9

## Installing required packages.
RUN apt-get update && apt-get install -y --force-yes \
    bedtools \
    libcurl4-gnutls-dev \
    libssh2-1-dev \
    libssl-dev \
    libxml2-dev \
    parallel \
    python \
    r-base \
    r-base-dev \
    wget \
    zlib1g-dev

## Make ATACseeker directory. 
RUN mkdir -p /atacseeker

## Copy scripts & reference to atacseeker folder
COPY atacseeker/scripts /atacseeker 
COPY atacseeker/reference /atacseeker

## Install R packages 
RUN Rscript /atacseeker/scripts/install_packages.R

## Install bedGraphToBigWig
RUN wget --directory-prefix=/tmp http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
RUN cp /tmp/bedGraphToBigWig /usr/local/bin

## Install RStudio for pandoc libraries, required for rmarkdown
## RStudio is removed once pandoc has been copied to bin
RUN wget --directory-prefix=/tmp https://download1.rstudio.org/rstudio-0.99.486-amd64-debian.tar.gz
RUN tar -zxvf /tmp/rstudio-0.99.486-amd64-debian.tar.gz -C /tmp && \
    cp /tmp/rstudio-0.99.486/bin/pandoc/* /bin && \
    rm -rf /tmp/rstudio*

## Install samtools v < 1.0
RUN wget --directory-prefix=/tmp https://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2
RUN tar -zxvf /tmp/samtools-0.1.19.tar.bz2 -C /tmp && \
    cd /tmp/samtools-0.1.19 && make && \
    cp /tmp/samtools-0.1.19/samtools /usr/local/bin && \
    rm -rf /tmp/samtools-0.1.19.tar.bz2
