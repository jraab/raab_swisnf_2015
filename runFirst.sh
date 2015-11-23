#!/bin/bash
# This script downloads some intermediate files that contain coverages for each of the ENCODE datasets over
# the peak sets used in the paper. 
# Running this will allow the figures to be generated using createFigs.sh
# It requires curl which should be on linux and osx systems 
 
mkdir -p data

# This pulls in the processed encode signal and enrichment over the ARID peaks
# Needed to create some of the metagene plots. 
curl -u swisnf:swisnf ftp://rc-ns-ftp.its.unc.edu/encode_coverages.tar.gz | tar -xvz -C output/

# This will bring in the data/ directory I used and includes most of my data and intermediate processed files
# I did not include fastq.gz, but the bam files that contain the aligned reads I used for counting and 
# visualization are included. 
# If you want fastq.gz raw reads please see GEO GSE69568
curl -u swisnf:swisnf ftp://rc-ns-ftp.its.unc.edu/data.tar.gz | tar -xvz -C data/

# To get the expression data from all encode cell lines run this line 
curdir=$(pwd)
cd data/external/encode_expr && xargs -n 1 curl -O -L < ../encode_expr.txta && cd $curdir


# Dependencies. 

#pip install deeptools

# Rscript code/utils/install_R_deps.R
