#!/usr/bin/sh
#$ -cwd 
#$ -o /magnuson-lab/jraab/analysis/swi_snf_final/code/logs/encodecoverage.err
#$ -e /magnuson-lab/jraab/analysis/swi_snf_final/code/logs/encodecoverage.err

#script is the same as the arid coverages but I'm using it to call a different set of bam files - txn factors and histones
module load samtools
module load python/2.7.6
source ~/virtualenvs/base/bin/activate 

echo $I
if [ "$I" == "NULL" ] ; then 
   echo a
   echo $OUT
   python ~/scripts/coverageFrame.py -p $P -b $B -binsize 10 -binaction mean -bamname $BNAME -pname $PNAME -o $OUT -m 

else 
   echo b
   python ~/scripts/coverageFrame.py -p $P -b $B -i $I -binsize 10 -binaction mean -bamname $BNAME -pname $PNAME -o $O -m 
fi
