#!/bin/bash
#$ -cwd 

module load bowtie2
module load samtools

PROJ='/magnuson-lab/jraab/analysis/swi_snf_final/'
OUTPUTDIR='/magnuson-lab/jraab/scratch/'
INPUTDIR='/magnuson-lab/jraab/analysis/swi_snf_final/data/chip/'
FINALDIR='/magnuson-lab/jraab/analysis/swi_snf_final/data/chip/processed/'

echo 'start'
if [ ! -d ${FINALDIR}merged ]; then
   mkdir -p ${FINALDIR}merged
fi 
echo $STUB

echo ${INPUTDIR}${STUB}
bowtie2 -p 8 --met-file $PROJ/code/logs/${STUB}_rep1.bowtie.stats -S ${OUTPUTDIR}${STUB}_rep1.sam -x /magnuson-lab/jraab/indexes/hg19 -U ${INPUTDIR}${STUB}_rep1.fastq.gz
bowtie2 -p 8 --met-file $PROJ/code/logs/${STUB}_rep2.bowtie.stats -S ${OUTPUTDIR}${STUB}_rep2.sam -x /magnuson-lab/jraab/indexes/hg19 -U ${INPUTDIR}${STUB}_rep2.fastq.gz 

convert all sam to bam 
FILES=$(find ${OUTPUTDIR} -name ${STUB}*.sam) 
echo $FILES
for IN in $FILES; do 
   OUT=$(basename ${IN} | awk -F. '{print $1}') 
   samtools view -bS $IN -o ${OUTPUTDIR}$OUT.bam
done

#sort bams
BAMS=$(find $OUTPUTDIR -name ${STUB}*.bam) 
echo $BAMS
for bam in $BAMS; do  
   OUTBAM=$(basename $bam| awk -F. '{print $1}')
   samtools sort $bam ${OUTPUTDIR}${OUTBAM}.sorted
done

#remove duplicates
BAMS=$(find $OUTPUTDIR -name ${STUB}*.sorted.bam)
for bam in $BAMS; do 
   OUTBAM=$(basename $bam |awk -F. '{print $1'})
   samtools rmdup $bam ${FINALDIR}${OUTBAM}.sorted.rmdup.bam
done

#merge bams
BAMS=$(find $FINALDIR -name ${STUB}*.sorted.rmdup.bam)
#echo $BAMS
OUTBAM=${OUTPUTDIR}${STUB}.rmdup.merged.bam
#samtools merge $OUTBAM $BAMS

#sort bam
echo $OUTBAM
samtools sort ${OUTPUTDIR}${STUB}.rmdup.merged.bam ${FINALDIR}${STUB}.rmdup.merged.sorted

#index bams 
samtools index ${FINALDIR}merged/${STUB}.rmdup.merged.sorted.bam

