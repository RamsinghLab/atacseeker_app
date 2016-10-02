#/bin/bash -xv

## freeBayes
cd /tmp
git clone --recursive git://github.com/ekg/freebayes.git
cd freebayes
make
make install

## vcflib
cd /tmp
git clone --recursive https://github.com/vcflib/vcflib.git
cd vcflib
make
cp bin/vcffilter /usr/local/bin

## bamaddrg
cd /tmp
git clone --recursive https://github.com/ekg/bamaddrg
cd bamaddrg
make
cp bamaddrg /usr/local/bin
