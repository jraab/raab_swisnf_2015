# Reviewers want a list of each gene occupied by each ARID and the consequence on trnascription. 
library(dplyr) 
library(tidyr) 
library(readr) 
source('code/misc/chippeak_annotation.R') # this contains some helper functions 

# Genes
arid1a_genes <- read.csv('output/diffexp/tables/arid1a_full_results.csv') %>% add_rownames()
arid1b_genes <- read.csv('output/diffexp/tables/arid1b_sig_results.csv') %>% add_rownames()
arid2_genes  <- read.csv('output/diffexp/tables/arid2_sig_results.csv') %>% add_rownames() 

# Peaks 
arid1a_peaks <- read.csv('output/macs_peaks/cleaned/arid1a_annotated.csv') 
arid1b_peaks <- read.csv('output/macs_peaks/cleaned/arid1b_annotated.csv') 
arid2_peaks  <- read.csv('output/macs_peaks/cleaned/arid2_annotated.csv') 

arid1a_peaks_to_genes <- arid1a_peaks %>% 
  group_by(nearestTSS) %>% 
  summarize(numberPeaks = n() ) %>%
  mutate(ARID = 'ARID1A') %>% 
  left_join(arid1a_genes, by = c('nearestTSS'='rowname')) %>%
  select(nearestTSS, numberPeaks, ARID, log2FoldChange, padj)
   
arid1b_peaks_to_genes <- arid1b_peaks %>% 
  group_by(nearestTSS) %>% 
  summarize(numberPeaks = n() ) %>%
  mutate(ARID = 'ARID1B') %>% 
  left_join(arid1a_genes, by = c('nearestTSS'='rowname')) %>%
  select(nearestTSS, numberPeaks, ARID, log2FoldChange, padj)

arid2_peaks_to_genes <- arid2_peaks %>%
  group_by(nearestTSS) %>% 
  summarize(numberPeaks = n() ) %>%
  mutate(ARID = 'ARID2') %>% 
  left_join(arid1a_genes, by = c('nearestTSS'='rowname')) %>%
  select(nearestTSS, numberPeaks, ARID, log2FoldChange, padj)

write_csv(arid1a_peaks_to_genes, 'output/diffexp/tables/arid1a_expchanges_directbound.csv') 
write_csv(arid1b_peaks_to_genes, 'output/diffexp/tables/arid1b_expchanges_directbound.csv') 
write_csv(arid2_peaks_to_genes,  'output/diffexp/tables/arid2_expchanges_directbound.csv') 

all <- rbind(arid1a_peaks_to_genes, arid1b_peaks_to_genes, arid2_peaks_to_genes)

g <- all %>% 
  ggplot(aes(x= log2FoldChange, y=  -log10(padj) ) ) + geom_point() + facet_wrap(~ARID)

ggsave('output/diffexp/plots/txn_changes_associated_withpeaks.pdf', g) 
