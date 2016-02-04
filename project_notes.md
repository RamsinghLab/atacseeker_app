# ATACseeker App #

## Mt-DNA Variant Calling ##

## ATACseq Analysis Pipeline ##

### preseqR ###

### csaw ###

### lola ###

## Docker ##

## Native App ##


**testing basespace app**

`ilmn` has a weird way of testing these apps. essentially, i downloaded the native app development [image](https://da1s119xsxmu0.cloudfront.net/sites/developer/native/nativeappsvm/BaseSpace%20Native%20Apps%20VM%20v14.ova) and then loaded it in virtual box. You need to login to the image via the terminal - `ssh basespace@localhost -p2222` . The password is `basespace`. 

This image contains all software to test native apps - docker, spacedock etc. Also the versions should be up-to-date with the ones `ilmn` has. 

Associated with every app, is a `spacedock` command which in my case is - `sudo spacedock -a d31c319955504b5ba162a51a0026d693 -m https://mission.basespace.illumina.com` . The idea is to run this command locally inside the VM and then trigger a test by doing a `Send To Local Agent`.
