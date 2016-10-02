#/bin/bash -xv

## bamaddrg
cd /tmp
git clone --recursive https://github.com/ekg/bamaddrg
cd bamaddrg
make
cp bamaddrg /usr/local/bin
cd /tmp
rm -rf bamaddrg

## freeBayes
cd /tmp
git clone --recursive git://github.com/ekg/freebayes.git
cd freebayes
make
make install
cd /tmp
rm -rf freebayes

## samblaster
cd /tmp
git clone https://github.com/GregoryFaust/samblaster
cd samblaster
make
cp samblaster /usr/local/bin
rm -rf samblaster

## vcflib
cd /tmp
git clone --recursive https://github.com/vcflib/vcflib.git
cd vcflib
make
cp bin/vcffilter /usr/local/bin
cd /tmp
rm -rf vcflib
