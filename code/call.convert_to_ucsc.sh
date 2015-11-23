#!bin/bash

#use this script to call a single qsub job for each mapping convesion -  when all jobs return run makehub.py outside a compute node 

PROJ=/magnuson-lab/jraab/analysis/swi_snf_final/
OUTRNA=${PROJ}data/rna/processed/
OUTCHIP=${PROJ}data/chip/processed/
LOGS=${PROJ}code/logs/


if [ ! -d ${OUTRNA}/bw ] ; then 
   mkdir -p ${OUTRNA}/bw
fi

if [ ! -d ${OUTCHIP}/bw ]; then 
   mkdir -p ${OUTCHIP}/bw
   mkdir -p ${OUTCHIP}/bb
fi 
<<<<<<< HEAD
CHIPSEQSTUB="input"
#jCHIPSEQSTUB="input arid1a arid1b arid2 snf5" 
=======

CHIPSEQSTUB="arid1a arid1b arid2 snf5" 
>>>>>>> e33a5d89c0f5d75afae01dde9452063e2287f8d6
#RNASEQSTUB="ntg arid1a arid1b arid2" 

for chip in $CHIPSEQSTUB; do 
   export STUB=$chip
   export OUTDIR=$OUTCHIP
   export TYPE='chip'
   export genome='/magnuson-lab/jraab/annotations/genome_hg19.txt'
   qsub -V  -cwd -o ${LOGS}${chip}.chip_convertucsc.out -e ${LOGS}${chip}.chip.convertucsc.err pipeline.convert_to_ucsc.sh
done


for rna in $RNASEQSTUB; do 
   export STUB=$rna
   export OUTDIR=$OUTRNA
   export TYPE='rna'
   qsub -V -cwd -o ${LOGS}${rna}.rna_convertucsc.out -e ${LOGS}${rna}.rna.convertucsc.err pipeline.convert_to_ucsc.sh
done
