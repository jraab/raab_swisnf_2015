# Define direct targets in specific categories. 
library(magrittr) 
library(dplyr)
library(ggplot2) 
source('code/theme_paper.R')
# load data for each arid bound set of peaks
arid1a_peaks <- read.csv('output/macs_peaks/cleaned/arid1a_annotated.csv')
arid1b_peaks <- read.csv('output/macs_peaks/cleaned/arid1b_annotated.csv') 
arid2_peaks <- read.csv('output/macs_peaks/cleaned/arid2_annotated.csv') 

# How many nearest genes expression changes in the RNAseq set
arid1a_sig <- read.csv('output/diffexp/tables/arid1a_sig_results.csv') %>%
  add_rownames() %>% select(rowname, log2FoldChange) 
arid1b_sig <- read.csv('output/diffexp/tables/arid1b_sig_results.csv') %>%
  add_rownames() %>% select(rowname, log2FoldChange)
arid2_sig <- read.csv('output/diffexp/tables/arid2_sig_results.csv') %>% 
  add_rownames() %>% select(rowname, log2FoldChange)

# Define direct targets that are significantly changed. 
arid1a_direct <- arid1a_sig[arid1a_sig$rowname %in% arid1a_peaks$nearestTSS, ]
arid1b_direct <- arid1b_sig[arid1b_sig$rowname %in% arid1b_peaks$nearestTSS, ]
arid2_direct  <- arid2_sig[arid2_sig$rowname %in% arid2_peaks$nearestTSS, ]

# Now define the comparison that are most interesting. 
# --------------------------------------------------------------------------

# Arid1a and Arid2 cooperative and competitive genes. 
combined <- inner_join(arid1a_direct, arid2_direct, by = 'rowname') %>% 
  mutate(type = ifelse(sign(log2FoldChange.x) == sign(log2FoldChange.y), 'coop', 'compete') )
print(table(combined$type)  )

# 14 competive interaction, 23 cooperative interactions 


# Arid1b and Arid2 cooperative and repressed 
arid1b_arid2 <- inner_join(arid1b_direct, arid2_direct, by = 'rowname') %>% 
  filter(sign(log2FoldChange.x) == sign(log2FoldChange.y) ) %>% 
  mutate(type = ifelse(log2FoldChange.x < 0, 'Activated', 'Repressed') )
table(arid1b_arid2$type ) 

# 60 % are repressed ( direct targets are lower than overall numbers) 
repressed <- arid1b_arid2 %>% filter(type == 'Repressed') %>% arrange(desc(log2FoldChange.x + log2FoldChange.y) )
# 30% of these genes are smad targets (35/114)
# repressed genes involved in celld evelopment, signal transduction , and apoptosis

activated <- arid1b_arid2 %>% filter(type == 'Activated') %>% arrange(log2FoldChange.x - log2FoldChange.y) 
# slight enrichment for metabolic processes and lipid functions
# genes down regulated in some cancers (prostate, pdac) and associated with liver specificity

# create a data frame of all genes X all enrichments - can use this to pick things that may be more specific for a 
# particular gene set. 
f <- read.csv('output/encode_coverages/enrichments/')