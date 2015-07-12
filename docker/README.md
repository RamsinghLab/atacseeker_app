** Contains Docker files for creating the ATACseq image **

To run the alignment, make sure to mount the "tmp" folder containing the data. 

Something like this should work:
*docker run -v ~/projects/atacseq_counts/docker/tmp:/atacseq_counts -t asifzubair/atacseq*


