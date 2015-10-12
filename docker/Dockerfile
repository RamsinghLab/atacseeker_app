# This Docker file builds an image for the alignment pipeline. 

# Might be better to use the base-image here: phusion/baseimage
# For now, using ubuntu:latest
# However, should use a particular version of the image for the sake of reproducibility.

FROM ubuntu
MAINTAINER Asif Zubair <asif.zubair@gmail.com>

ENV DEBIAN_FRONTEND noninteractive

# Getting ready to install R!
RUN echo 'deb http://cran.stat.ucla.edu//bin/linux/ubuntu trusty/' | tee -a /etc/apt/sources.list.d/r.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9

# Installing required packages.
RUN apt-get update && apt-get install -y \
    bwa \
    git \
    r-base \
    r-base-dev \
    samtools

# Make ATACseeker directories. 
RUN mkdir -p /atacseq_counts/data /atacseq_counts/reference /atacseq_counts/scripts /atacseq_counts/install

COPY tmp/data /atacseq_counts/data
COPY tmp/reference /atacseq_counts/reference
COPY tmp/scripts /atacseq_counts/scripts

#COPY tmp/scripts tmp

# Run samblaster installation and run alignment script. 
RUN bash /atacseq_counts/scripts/install_samblaster.sh

#RUN bash /tmp/install_samblaster.sh 
CMD ["bash","/atacseq_counts/scripts/my_alignment_pipeline.sh"]


