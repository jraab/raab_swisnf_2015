#!/bin/bash

#script used to submit a pipeline job for each group of chipseq samples

#chipseq
CHIP="arid1a arid1b arid2 snf5" 
ERRORPATH=/magnuson-lab/jraab/analysis/swi_snf_final/code/logs/
for chip in $CHIP; do 
   echo $chip
   export STUB=$chip
   qsub -V -m abe -M jesse.r.raab@gmail.com -e ${ERRORPATH}${chip}.err \
   -o ${ERORRPATH}${chip}.log pipeline.chipreadhandler.sh
done

<<<<<<< HEAD
qsub -V -m abe -M jesse.r.raab@gmail.com -e ${ERRORPATH}.input.err \
   -o ${ERRORPATH}$input.log pipeline.input.sh

=======
>>>>>>> e33a5d89c0f5d75afae01dde9452063e2287f8d6




