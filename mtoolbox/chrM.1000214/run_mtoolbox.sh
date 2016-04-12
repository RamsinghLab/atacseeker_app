#!/bin/bash

#MToolBox.sh -i fastq -r RCRS

HOME="/home/cmb-07/sn1/asifzuba"

MToolBox.sh -i bam -r RCRS \
-m "-g ${HOME}/software_frozen/local/bin/gsnap -D ${HOME}/software_frozen/local/share/gmapdb/ -M chrM -H hg19RCRS " \
-a "-r ${HOME}/software_frozen/local/share/genomes/ -f chrM.fa -a hg19RCRS.fa -s ${HOME}/software_frozen/src/anaconda/envs/mtoolbox/bin/samtools " \
-c "-m ${HOME}/software_frozen/src/anaconda/envs/mtoolbox/bin/muscle" \
-s ${HOME}/software_frozen/src/anaconda/envs/mtoolbox/bin/samtools
