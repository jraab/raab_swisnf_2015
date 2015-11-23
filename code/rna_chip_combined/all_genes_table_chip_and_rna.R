# Script to generate table of ARID binding and gene expression change and adjusted pvalue for 
# each ARID

library(dplyr) 
library(tidyr) 
library(readr) 
library(gridExtra) 
library(ggplot2) 
source('code/misc/chippeak_annotation.R') # this contains some helper functions 
source('code/util/theme_paper.R')
source('code/function_defs/rnaseq_fig_fun.R')

# All genes 
all_genes <- read.csv('output/diffexp/tables/arid1a_full_results.csv') %>% add_rownames() %>%
  select(rowname) %>% tbl_df()

# Genes
arid1a_genes <- read.csv('output/diffexp/tables/arid1a_full_results.csv') %>% add_rownames()%>%
  select(rowname, log2FoldChange, padj) 
arid1b_genes <- read.csv('output/diffexp/tables/arid1b_full_results.csv') %>% add_rownames()%>%
  select(rowname, log2FoldChange, padj) 
arid2_genes  <- read.csv('output/diffexp/tables/arid2_full_results.csv') %>% add_rownames() %>% 
  select(rowname, log2FoldChange, padj)  

# Peaks 
arid1a_peaks <- read.csv('output/macs_peaks/cleaned/arid1a_annotated.csv') %>% add_rownames() 
arid1b_peaks <- read.csv('output/macs_peaks/cleaned/arid1b_annotated.csv') %>% add_rownames() 
arid2_peaks  <- read.csv('output/macs_peaks/cleaned/arid2_annotated.csv') %>% add_rownames()

# For each genes add column if it is in the peak list. 
all_genes$Arid1a_peak <- ifelse(all_genes$rowname %in% arid1a_peaks$nearestTSS, 1, 0) 
all_genes$Arid1b_peak <- ifelse(all_genes$rowname %in% arid1b_peaks$nearestTSS, 1, 0)
all_genes$Arid2_peak <- ifelse(all_genes$rowname %in% arid2_peaks$nearestTSS, 1, 0) 

# Merge data from each peak set 
all_genes <- all_genes %>% left_join(arid1a_genes, by = 'rowname') %>% 
  left_join(arid1b_genes, by = 'rowname') %>%
  left_join(arid2_genes, by = 'rowname') 
colnames(all_genes) <- c('Gene', 'ARID1A_peak', 'ARID1B_peak', 'ARID2_peak', 'ARID1A_Log2FC', 'ARID1A_padj', 
                         'ARID1B_Log2FC', 'ARID1B_padj', 'ARID2_Log2FC', 'ARID2_padj')

# _
write_csv(all_genes, 'output/supplemental_data/all_genes_chip_and_expression.csv')
