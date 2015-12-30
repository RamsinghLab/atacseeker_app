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
##    bwa \
##    git \
    libcurl4-gnutls-dev \
    libxml2-dev \
    r-base \
    r-base-dev \
##    samtools \
    wget

## Make ATACseeker directories. 

## RUN mkdir -p /atacseeker/reference /atacseeker/scripts 
RUN mkdir -p /atacseeker/scripts

## Copy reference and scripts to atacseeker folder

## COPY tmp/reference /atacseeker/reference
COPY tmp/scripts /atacseeker/scripts

## Install R packages 

RUN Rscript /atacseeker/scripts/install_packages.R

## Install RStudio for pandoc libraries
## Required for rmarkdown
## RStudio is removed once pandoc has been copied to bin

RUN wget --directory-prefix=/tmp https://download1.rstudio.org/rstudio-0.99.486-amd64-debian.tar.gz
RUN tar -zxvf /tmp/rstudio-0.99.486-amd64-debian.tar.gz -C /tmp
RUN cp /tmp/rstudio-0.99.486/bin/pandoc/* /bin
RUN rm -rf /tmp/rstudio*

## Install samblaster from Github repository 

## RUN bash /atacseeker/scripts/install_samblaster.sh

## Build the index

## RUN bwa index -a bwtsw /atacseeker/reference/hg19.fa.gz

## Render the Rmarkdown file

## CMD Rscript -e 'rmarkdown::render("/atacseeker/scripts/atacseeker.Rmd")'