#!/bin/bash

# This script will download samblaster and install it
# It then adds copies samblaster to the local/bin so that it can be called easily

cd /Package_content/software
git clone git://github.com/GregoryFaust/samblaster.git

cd samblaster
make

cp samblaster /usr/local/bin/.
