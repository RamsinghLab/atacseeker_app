##
## This Docker file builds an image for the ATACseeker pipeline. 
##

FROM ubuntu:14.04
MAINTAINER Asif Zubair <asif.zubair@gmail.com>

ENV DEBIAN_FRONTEND noninteractive

# Getting ready to install R

RUN echo 'deb http://cran.stat.ucla.edu//bin/linux/ubuntu trusty/' | tee -a /etc/apt/sources.list.d/r.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9

# Installing required packages.

RUN apt-get update && apt-get install -y \
    bwa \
    git \
    r-base \
    samtools \
    wget

# Make ATACseeker directories. 

RUN mkdir -p /atacseker/reference /atacseeker/scripts /atacseeker/install

# Copy reference and scripts to atacseeker folder

COPY tmp/reference /atacseeker/reference
COPY tmp/scripts /atacseeker/scripts

# Install R packages -  knitr, rmarkdown, statmod, locfit, edgeR, csaw

RUN R CMD BATCH /atacseeker/scripts/install_packages.R

# Install RStudio for pandoc libraries
# Required for rmarkdown
# RStudio is removed once pandoc has been copied to bin

RUN wget -q https://download1.rstudio.org/rstudio-0.99.486-amd64-debian.tar.gz /tmp
RUN tar -zxf /tmp/rstudio-0.99.486-amd64-debian.tar.gz
RUN cp /tmp/rstudio-0.99.486/bin/pandoc/* /bin
RUN rm -rf /tmp/rstudio*

# Install samblaster from Github repository 

RUN bash /atacseeker/scripts/install_samblaster.sh

# Render the Rmarkdown file

CMD Rscript -e 'rmarkdown::render("/atacseeker/scripts/atacseeker.Rmd")'
