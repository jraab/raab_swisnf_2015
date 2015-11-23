#!/bin/bash 

#script to convert all chipseq and rnaseq files to ucsc .bw tracks 
#also converts peak files to bb for viewing on ucsc 
# last part will convert the bw and bb fiels from all this to a track hub and upload to kure for viewing. 
PROJ=/magnuson-lab/jraab/analysis/swi_snf_final/
SCRATCH=$HOME/scratch/

module load python/2.7.6
module load samtools
module load bedtools 
PATH=$PATH:/magnuson-lab/jraab/bin
source /magnuson-lab/jraab/virtualenvs/base/bin/activate 
mkdir -p ${SCRATCH}/rna
mkdir -p ${SCRATCH}/chip

echo $STUB

if [ $TYPE == 'chip' ] ; then   
   if [ $STUB == 'input' ]; 
   then 
      bamfile=$(find ${OUTDIR}/merged -iname ${STUB}.sorted.bam)
      python ${HOME}/scripts/bamToWig_chip.py -b $bamfile -g hg19 --norm -d ${SCRATCH}/chip/ -frag 400
   else 
      bamfile=$(find ${OUTDIR}/merged -iname ${STUB}*.rmdup.merged.sorted.bam)   
      python $HOME/scripts/bamToWig_chip.py -b $bamfile -g hg19 --norm -d $SCRATCH/chip/ -frag 400
   fi 
   echo $bamfile
   bgfile=$(find ${HOME}/scratch/chip/ -name ${STUB}.bg ) 
   echo $bgfile
   NAME=$(basename $bgfile | awk -F. '{print $1}')  
   sort -k1,1 -k2,2n $bgfile > ${SCRATCH}chip/${NAME}.sorted.tmp
   head -n -1 ${SCRATCH}chip/${NAME}.sorted.tmp > ${SCRATCH}chip/${NAME}.sorted.tmp2
   bedGraphToBigWig ${SCRATCH}chip/${NAME}.sorted.tmp2 $genome $OUTDIR/bw/${NAME}.bw 
   
   
   if [ $STUB != 'input' ]; then 
      bedfile=/magnuson-lab/jraab/analysis/swi_snf_final/output/macs_peaks/cleaned/${STUB}_peaks_cf.bed
      sort -k1,1 -k2,2n $bedfile > ${SCRATCH}chip/${NAME}.peaks.sorted.tmp
      bedToBigBed ${SCRATCH}chip/${NAME}.peaks.sorted.tmp $genome ${OUTDIR}/bb/${NAME}.bb
   fi
   #rm ${SCRATCH}chip/${NAME}*tmp*
               
fi 

if [ $TYPE == 'rna' ]; then  
   bamfile=$(find ${OUTDIR}/merged -iname ${STUB}*.merged.sorted.bam)  
   echo $OUTDIR
   bamfile=$(find ${OUTDIR}/merged -iname ${STUB}*.rmdup.merged.sorted.bam) 
   echo $bamfile
   python $HOME/scripts/bamToWig_chip.py -b $bamfile -g hg19 --norm -d $SCRATCH/chip/ -frag 250
   bgfile=$(find ${HOME}/scratch/chip/ -name ${STUB}.bg ) 
   echo $bgfile
   NAME=$(basename $bgfile | awk -F. '{print $1}')
   echo ${SCRATCH}chip/${NAME}
   sort -k1,1 -k2,2n $bgfile > ${SCRATCH}chip/${NAME}.sorted.tmp
   head -n -1 ${SCRATCH}chip/${NAME}.sorted.tmp > ${SCRATCH}chip/${NAME}.sorted.tmp2
   bedGraphToBigWig ${SCRATCH}chip/${NAME}.sorted.tmp2 $genome $OUTDIR/bw/${NAME}.bw 
   bedfile=/magnuson-lab/jraab/analysis/swi_snf_final/output/macs_peaks/cleaned/${STUB}_peaks_cf.bed
   sort -k1,1 -k2,2n $bedfile > ${SCRATCH}chip/${NAME}.peaks.sorted.tmp
   bedToBigBed ${SCRATCH}chip/${NAME}.peaks.sorted.tmp $genome ${OUTDIR}/bb/${NAME}.bb
   rm ${SCRATCH}chip/${NAME}*tmp*
               
fi 



