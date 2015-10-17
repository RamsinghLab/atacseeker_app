#!/bin/bash

####
## This script will download samblaster and install it
## It then copies samblaster to the local/bin so that it can be called easily
####

cd /tmp
git clone git://github.com/GregoryFaust/samblaster.git

cd samblaster
make

cp samblaster /usr/local/bin/
