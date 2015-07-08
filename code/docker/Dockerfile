# This is Docker file builds an image for the alignment pipeline. 

# Might be better to use the base-image here: phusion/baseimage
# For now, using ubuntu:latest
# However, best practices inidicate to use a particular version of the image for the sake of reproducibility.

FROM ubuntu
MAINTAINER Asif Zubair <asif.zubair@gmail.com>

# Honestly, I have no idea why i am running this ENV command, must ask Anthony.
ENV DEBIAN_FRONTEND noninteractive

# Getting ready to install R!
RUN echo 'deb http://cran.stat.ucla.edu//bin/linux/ubuntu trusty/' | sudo tee -a /etc/apt/sources.list.d/r.list
RUN sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9

# Best practices indicate not to run apt-get update on a single line and have multi line installation
RUN apt-get update && install -y \
    r-base \
    bwa \
    samtools \
    git

#COPY tmp/software /bin
COPY tmp/samplefastq /Package_data/sample_fastq

CMD [./my_alignment_script.sh]
