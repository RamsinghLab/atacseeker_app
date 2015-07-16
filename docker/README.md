**Contains Docker files for creating the ATACseq image**

Running the command: 

`docker build -t image-name .`

in this folder will create a docker image with the following software installed:

bwa, git, R, samblaster, samtools


To run the alignment inside the docker container, directly run the pipeline by doing this:

`docker run -t image-name`


Previously, was trying to not copy data & reference to image. In this case, make sure to mount the "tmp" folder containing the data. 

Something like this should work:

`docker run -v ~/projects/atacseq_counts/docker/tmp:/atacseq_counts -t image-name`
