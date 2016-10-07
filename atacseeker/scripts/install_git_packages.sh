#/bin/bash -xv

TMP_DIR="/tmp"

## bamaddrg
cd ${TMP_DIR}
git clone --recursive https://github.com/ekg/bamaddrg
cd bamaddrg
make
cp bamaddrg /usr/local/bin
cd ${TMP_DIR}
rm -rf bamaddrg

## freeBayes
cd ${TMP_DIR}
git clone --recursive git://github.com/ekg/freebayes.git
cd freebayes
make
make install
cd ${TMP_DIR}
rm -rf freebayes

## samblaster
cd ${TMP_DIR}
git clone https://github.com/GregoryFaust/samblaster
cd samblaster
make
cp samblaster /usr/local/bin
cd ${TMP_DIR}
rm -rf samblaster

## vcflib
cd ${TMP_DIR}
git clone --recursive https://github.com/vcflib/vcflib.git
cd vcflib
make
cp bin/vcffilter /usr/local/bin
cd ${TMP_DIR}
rm -rf vcflib
