#!/bin/bash
#$ -cwd 
#$ -e /magnuson-lab/jraab/analysis/swi_snf_final/code/logs/fixpeaks.err
#$ -o /magnuson-lab/jraab/analysis/swi_snf_final/code/logs/fixpeaks.out

module load python/2.7.6
module load bedtools 
source /magnuson-lab/jraab/virtualenvs/base/bin/activate 
OUTDIR=/magnuson-lab/jraab/analysis/swi_snf_final/output/macs_peaks/cleaned/
if [ ! -d ${OUTDIR} ] ; then 
   mkdir -p ${OUTDIR} 
fi 
python /magnuson-lab/jraab/analysis/swi_snf_final/code/chip_analysis/merge_and_filter_peaks.py
