#!/bin/sh
#$ -cwd 
#$ -o ../logs/k27ac_peaks.out
#$ -e ../logs/k27ac_peaks.err

module load python/2.7.6
module load bedtools
source ~/virtualenvs/base/bin/activate 

K27=/magnuson-lab/jraab/ENCODE/datafiles/HepG2/combined/wgEncodeBroadHistoneHepg2H3k27acStdAln.sorted.merged.bam
INPUT=/magnuson-lab/jraab/ENCODE/datafiles/HepG2/combined/wgEncodeBroadHistoneHepg2ControlStdAln.sorted.merged.bam
OUT=/magnuson-lab/jraab/ENCODE/peaks/
SWI=/magnuson-lab/jraab/analysis/swi_snf_final/output/
if [[ ! -d $OUT ]]; then 
   mkdir -p $OUT
fi

# call peaks
macs2 callpeak -t $K27 -c $INPUT -n $OUT/k27ac_hepg2 -g 2.7e9 

# get coverage in each peaks - easier to merge these in R as desired
bedtools coverage -abam $K27 -b ${OUT}/k27ac_hepg2_peaks.narrowPeak > ${SWI}k27ac_coverage.btout.txt
