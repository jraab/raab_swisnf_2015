#!/bin/bash

#script used to submit a pipeline job for each group of chipseq samplse
#CHIP="arid1a arid1b arid2 snf5" 
RNA="ARID1A ARID1B ARID2 NTG"
ERRORPATH=/magnuson-lab/jraab/analysis/swi_snf_final/code/logs/rna/
if [ ! -d $ERRORPATH ] ; then 
   mkdir -p $ERRORPATH
fi 

#rnasq
for rna in $RNA; do 
   export STUB=$rna
   qsub -V -m abe -M jesse.r.raab@gmail.com -e ${ERRORPATH}${rna}.rna.err \
   -o ${ERRORPATH}${rna}.rna.log pipeline.rnareadhandler.sh
done

