#!/bin/bash
#$ -cwd
#$ -e /magnuson-lab/jraab/analysis/swi_snf_final/code/logs/chipseq_analyis_pipeline.err
#$ -o /magnuson-lab/jraab/analysis/swi_snf_final/code/logs/chipseq_analysis_pipeline.err

# Notes
# ----------------------------------------------------------
# Script shoudl be run directly  as
# qsub chipseq_analysis_pipeline.sh
# scripts are fast, so no need to pipe out each one idividually

# load modules
# ----------------------------------------------------------
module load r
module load python/2.7.6
source /magnuson-lab/jraab/virtualenvs/base/bin/activate

<<<<<<< HEAD
# ----------------------------------------------------------
# Annotate Arid peak with genomic bound regions 
Rscript --vanilla chipseq_analysis/annotate_arid_peaks

#-----------------------------------------------------------
# Arid overlaps (against snf5 and each other 
Rscript --vanilla chipseq_analysis/arid_overlap.R


#----------------------------------------------------------
# Plot histone modification enrichment at Arids

#----------------------------------------------------------
# Plot transcription factor modifications at Arids


=======

# ----------------------------------------------------------
# Script to take macs chipseq peaks and convert to final data

# ----------------------------------------------------------
# Clean macs called peaks - remove blacklist, merge nearby



# ----------------------------------------------------------
# Annotate Arid peask with genomic bound regions 
>>>>>>> e33a5d89c0f5d75afae01dde9452063e2287f8d6
