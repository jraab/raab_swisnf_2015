#!/bin/bash 

FILES=$(find data/rna/processed/merged -name "*.merged.sorted.bam")
echo $FILES
BED=$(find output/misc/ -name "*_enhancers.bed")
echo $BED
OUTPUT='output/rnacov/'
if [[ ! -d $OUTPUT ]]; then 
   mkdir -p $OUTPUT
fi

for i in $FILES; do 
   for j in $BED; do 
      F=$(basename $i | awk -F. '{print $1}')
      J=$(basename $j | awk -F. '{print $1}')
      bedtools coverage -abam $i -b $j -counts > ${OUTPUT}${F}_${J}.txt 
      
   done
done
