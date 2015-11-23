#!/bin/bash
#$ -cwd 

module load tophat
module load bowtie2
module load samtools
module load python/2.7.6
source /magnuson-lab/jraab/virtualenvs/base/bin/activate

PROJ=/magnuson-lab/jraab/analysis/swi_snf_final/
OUTPUTDIR=/magnuson-lab/jraab/scratch/rna/
INPUTDIR=/magnuson-lab/jraab/analysis/swi_snf_final/data/rna/
FINALDIR=/magnuson-lab/jraab/analysis/swi_snf_final/data/rna/processed/
ANNO=/magnuson-lab/shared/jraab/annotations/gencode.v16.annotation.gtf

if [ ! -d $FINALDIR ]; then
   mkdir -p $FINALDIR
fi 

echo $STUB
ALLFILES=$(find $INPUTDIR -name $STUB*.fastq.gz) 

for f in $ALLFILES; do 
   #get the rep number
    REP=$(basename $f | awk -F_ '{print $2}' | awk -F. '{print $1}' )
   echo ${INPUTDIR}${STUB}_${REP}
   echo $f
   tophat -p 8 --library-type fr-firststrand --output-dir ${OUTPUTDIR}${STUB}_${REP}/ /magnuson-lab/jraab/indexes/hg19 $f
   #sort bams
   bam=$(find ${OUTPUTDIR}${STUB}_${REP}/ -name accepted_hits.bam)   
   samtools sort $bam ${FINALDIR}${STUB}_${REP}.sorted
   done


if [ ! -d ${FINALDIR}/merged ] ; then 
   mkdir -p ${FINALDIR}/merged 
fi 

#merge bams
SORTEDBAMS=$(find $FINALDIR -name ${STUB}*.sorted.bam)
echo $SORTEDBAMS
echo $STUB
OUTBAM=${OUTPUDIR}${STUB}.merged.bam
samtools merge -h ${FINALDIR}${STUB}_1.sorted.bam ${OUTBAM} $SORTEDBAMS

#sort merged bams
BAMS=$(find $FINALDIR -name ${STUB}*.merged.bam) 
for bam in $BAMS; do 
   OUTBAM=$(basename $bam | awk -F. '{print $1}') 
   samtools sort $bam ${FINALDIR}${OUTBAM}.merged.sorted
done

#index bams 
for bam in $(find $FINALDIR -name ${STUB}*.sorted.bam) ; do
   samtools index $bam
done


if [ ! -d $FINALDIR/counts ] ; then
   mkdir -p $FINALDIR/counts
fi

#count reads 
BAMS=$(find ${FINALDIR}${STUB}_*.sorted.bam) 
echo $BAMS
for bam in ${BAMS};  do 
   NAME=$(basename $bam | awk -F. '{print $1}' )
   echo $NAME
   OUT=${FINALDIR}/counts/${NAME}.counts.tab
   samtools view ${bam} | htseq-count -s reverse -m intersection-strict -i gene_name -  ${ANNO} > ${OUT}
done
