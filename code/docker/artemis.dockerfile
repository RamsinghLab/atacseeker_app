#this is a dockerfile to run the command Artemis containerization

FROM ubuntu
MAINTAINER anthonycolombo60@gmail.com
ENV DEBIAN_FRONTEND noninteractive
RUN echo 'deb http://cran.stat.ucla.edu//bin/linux/ubuntu trusty/' | sudo tee -a /etc/apt/sources.list.d/r.list
RUN sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
RUN apt-get update 
RUN apt-get install r-base -y
COPY software /bin
COPY transcriptomes /Package_data/transcriptomes
COPY Samplefastq /Package_data/sample_fastq
