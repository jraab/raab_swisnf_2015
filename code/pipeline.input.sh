#!/bin/sh
#$ -cwd 

module load bowtie2
module load samtools

PROJ='/magnuson-lab/jraab/analysis/swi_snf_final/'
OUTPUTDIR='/magnuson-lab/jraab/scratch/'
INPUTDIR='/magnuson-lab/jraab/analysis/swi_snf_final/data/chip/'
<<<<<<< HEAD
FINALDIR='/magnuson-lab/jraab/analysis/swi_snf_final/data/chip/processed/merged'
=======
FINALDIR='/magnuson-lab/jraab/analysis/swi_snf_final/data/processed/'
>>>>>>> e33a5d89c0f5d75afae01dde9452063e2287f8d6

if [ ! -d $FINALDIR ]; then
   mkdir -p $FINALDIR
fi 
echo $STUB
STUB='input'
bowtie2 -p 8 --met-file $PROJ/code/logs/${STUB}.bowtie.stats -S ${OUTPUTDIR}${STUB}_rep0.sam -x /magnuson-lab/jraab/indexes/hg19 -U ${INPUTDIR}${STUB}_rep0.fastq.gz

#convert all sam to bam 
FILES=$(find ${OUTPUTDIR} -name ${STUB}*.sam) 
echo $FILES
for IN in $FILES; do 
   samtools view -bS $IN -o ${OUTPUTDIR}${STUB}.bam
done

#sort bams
BAMS=$(find $OUTPUTDIR -name ${STUB}*.bam) 
echo $BAMS
for bam in $BAMS; do  
   OUTBAM=$(basename $bam| awk -F. '{print $1}')
   samtools sort $bam ${FINALDIR}${OUTBAM}.sorted
done


#index bams 
for bam in $(find $FINALDIR -name ${STUB}*.sorted.bam) ; do
   samtools index $bam
done



