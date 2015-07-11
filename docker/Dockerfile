# This Docker file builds an image for the alignment pipeline. 

# Might be better to use the base-image here: phusion/baseimage
# For now, using ubuntu:latest
# However, best practices inidicate to use a particular version of the image for the sake of reproducibility.

FROM ubuntu
MAINTAINER Asif Zubair <asif.zubair@gmail.com>

ENV DEBIAN_FRONTEND noninteractive

# Getting ready to install R!
RUN echo 'deb http://cran.stat.ucla.edu//bin/linux/ubuntu trusty/' | sudo tee -a /etc/apt/sources.list.d/r.list
RUN sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9

# Best practices indicate not to run apt-get update on a single line and to have multi-line installation sorted alphabetically.
RUN apt-get update && apt-get install -y \
    bwa \
    git \
    r-base \
    samtools

# Make project directories. 
RUN mkdir -p /atacseq_counts/data /atacseq_counts/reference /atacseq_counts/scripts

#COPY tmp/data /atacseq_counts/data
#COPY tmp/reference /atacseq_counts/reference
COPY tmp/scripts /atacseq_counts/scripts

RUN bash /atacseq_counts/scripts/install_sambalster.sh
#CMD [./my_alignment_script.sh]
