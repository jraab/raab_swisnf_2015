#!/bin/bash
# create a bed for use with deeptools where each arid peak file is 
# filtered for the alone peaks
# and the mulitple peaks are concatenated at the end 
# should be indicated as 
#  #groupname 
# valid bed


arid1a=../../output/macs_peaks/cleaned/arid1a_annotated.csv
arid1b=../../output/macs_peaks/cleaned/arid1b_annotated.csv
arid2=../../output/macs_peaks/cleaned/arid2_annotated.csv

names=('Arid1a' 'Arid1b' 'Arid2')

k=0

out=../../output/macs_peaks/cleaned/all_peaks_fordt.bed
rm $out
for i in $arid1a $arid1b $arid2; do  
   
   cat $i | grep 'Single'| cut -d, -f 1,2,3,4 | tr "," "\t" | \
      sed 's/"//g' >> $out
   echo "#${names[k]}" >> $out
   k=$k+1
done

cat ../../output/macs_peaks/cleaned/multi_bound_peaks.csv | \
   cut -d, -f 1,2,3,4 | tr "," "\t" | \
   sed 's/"//g' >> $out
   echo "#Multiple" >> $out
