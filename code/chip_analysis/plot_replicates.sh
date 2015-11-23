#!/bin/bash
module load r
stubs="arid1a arid1b arid2 snf5" 
for s in $stubs; do 
   Rscript --vanilla plot_replicate_comparisons.R $s
done
