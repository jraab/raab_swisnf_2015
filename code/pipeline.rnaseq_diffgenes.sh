#!/bin/bash
#$ -cwd
#$ -e /magnuson-lab/jraab/analysis/swi_snf_final/code/logs/diff_expr.err
#$ -o /magnuson-lab/jraab/analysis/swi_snf_final/code/logs/diff_expr.out
module load python/2.7.6
module load r
source /magnuson-lab/jraab/virtualenvs/base/bin/activate 

COUNTDIR=/magnuson-lab/jraab/analysis/swi_snf_final/data/rna/processed/counts/

#step 1 - merge all count data from searpate RNA seq into 1 tabel 
#input = all files in COUNTDIR and a python script to process
#output = table of expression changes
/usr/bin/env python /magnuson-lab/jraab/analysis/swi_snf_final/code/rna_analysis/process_htseq_count.py -d ${COUNTDIR}

#step 2 - use R and Deseq2 to id diff expr genes 
#input table of count data from step 1 ; R script to process
#output table for each ARID giving all info for all genes 
#       table of normalized RPKM type value for each gene (rows) for each arid ( column ) 
#OUTPUTDIR=/magnuson-lab/jraab/analysis/swi_snf_final/output/diffexpr/
<<<<<<< HEAD
Rscript --vanilla differential_expression_test.R 
Rscript --vanilla rnaseq_figures.R

=======
#Rscript --vanilla diff_expr_arids ${OUTPUTDIR}
>>>>>>> e33a5d89c0f5d75afae01dde9452063e2287f8d6

# step 3 generate plots for figure 1 
# figures will be dumped in  $PROJ/output/diffexp/tables or 
#Rscript --vanilla rnaseq_figures.R
