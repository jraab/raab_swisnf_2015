#!/bin/bash
module load bedtools 
STUBS="arid1a arid1b arid2 snf5"
tmp=/magnuson-lab/jraab/scratch/
output=/magnuson-lab/jraab/analysis/swi_snf_final/data/chip/processed/replicate_coverages/
if [ ! -d $output ] ; then 
   mkdir -p $output
fi

if [ ! -e ${tmp}/hg19.windows.10kb.bed ] ; then 
   bedtools makewindows -g /magnuson-lab/jraab/annotations/genome_hg19.txt -w 10000 > $tmp/hg19.windows.10kb.bed
fi 

fdir=/magnuson-lab/jraab/analysis/swi_snf_final/data/chip/processed/
for stub in $STUBS; do 
   FILES=$(find $fdir -name ${stub}*rep*.sorted.bam)
   for f in $FILES; do 
      echo $f
      outname=$(basename $f | awk -F. '{ print $1}')
      printf "#!/bin/bash
              #$ -cwd        
              module load bedtools
              bedtools coverage -abam ${f} -b ${tmp}/hg19.windows.10kb.bed" > tmp.sh      
      cat tmp.sh
      qsub -e ${outname}.err -o ${output}${outname}.10kbcov.txt -N ${outname}.cov tmp.sh
      rm tmp.sh
   done
done


