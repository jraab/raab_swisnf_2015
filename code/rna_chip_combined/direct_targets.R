# expression values at genes bound by each ARID at promoter (NTG vs siARID) 
library(magrittr) 
library(dplyr)
library(ggplot2) 
source('code/util/theme_paper.R')
# Define expression values for each condition 
#-------------------------------------------------------------------------------------
countdata <- read.csv('output/diffexp/tables/count_data_rpkm.csv')
countdata %<>%
  group_by(Sample, gene) %>% 
  summarise(avg_value = mean(Value, na.rm =T), 
            sd = sd(Value, na.rm =T) ) 

#make a separate frame for each comparison
arid1a_comp <- countdata[countdata$Sample %in% c('ARID1A', 'NTG'), ]
arid1b_comp <- countdata[countdata$Sample %in% c('ARID1B', 'NTG'), ]
arid2_comp <- countdata[countdata$Sample %in% c('ARID4', 'NTG'), ]

# load data for each arid bound set of peaks
arid1a_peaks <- read.csv('output/macs_peaks/cleaned/arid1a_annotated.csv')
arid1b_peaks <- read.csv('output/macs_peaks/cleaned/arid1b_annotated.csv') 
arid2_peaks <- read.csv('output/macs_peaks/cleaned/arid2_annotated.csv') 

# create data frame of gene expression info for each arid with genes that are either
# bound or not by the arid

arid1a_bound <- arid1a_comp[arid1a_comp$gene %in% arid1a_peaks$nearestTSS, ] %>% mutate(type='bound')
arid1a_unbound <- arid1a_comp[!arid1a_comp$gene %in% arid1a_peaks$nearestTSS, ] %>% mutate(type = 'unbound') 
c_arid1a <- rbind(arid1a_bound, arid1a_unbound) 
c_arid1a %>% 
  ggplot(aes(x= Sample, y = log(avg_value), fill = Sample, color = Sample) ) + geom_jitter(alpha  = 0.2) + facet_wrap(~type, scales = 'free')

# other way to look is to plot log2 fold change for the different sets of genes against one another

arid1a_all <- read.csv('output/diffexp/tables/arid1a_full_results.csv') %>% 
  add_rownames() %>% select(rowname, log2FoldChange) %>% mutate(arid = 'Arid1a')
arid1b_all <- read.csv('output/diffexp/tables/arid1b_full_results.csv') %>% 
  add_rownames() %>% select(rowname, log2FoldChange) %>% mutate(arid = 'Arid1b')
arid2_all <- read.csv('output/diffexp/tables/arid2_full_results.csv') %>% 
  add_rownames() %>% select(rowname, log2FoldChange) %>% mutate(arid = 'Arid2')

all <- rbind(arid1a_all, arid1b_all, arid2_all) 

arid1a_bound <- all[all$rowname %in% arid1a_peaks$nearestTSS, ]
arid1a_bound %>% ggplot(aes(arid,log2FoldChange)) +geom_violin()

selcrit <- arid2_peaks$nearestTSS %>% filter(as.character(type) == 'Promoter') 
arid2_bound <- all[all$rowname %in% selcrit, ]
arid2_bound %>% ggplot(aes(x = arid, y = log2FoldChange)) + geom_density()


# How many nearest genes expression changes in the RNAseq set
arid1a_sig <- read.csv('output/diffexp/tables/arid1a_sig_results.csv') %>%
  add_rownames() %>% select(rowname, log2FoldChange) 
arid1b_sig <- read.csv('output/diffexp/tables/arid1b_sig_results.csv') %>%
  add_rownames() %>% select(rowname, log2FoldChange)
arid2_sig <- read.csv('output/diffexp/tables/arid2_sig_results.csv') %>% 
  add_rownames() %>% select(rowname, log2FoldChange)

arid1a_sig %<>% mutate(direct = rowname %in% arid1a_peaks$nearestTSS, arid = 'Arid1a') 
arid1b_sig %<>% mutate(direct = rowname %in% arid1b_peaks$nearestTSS, arid = 'Arid1b') 
arid2_sig %<>% mutate(direct = rowname %in% arid2_peaks$nearestTSS, arid = 'Arid2') 
all_sig <- rbind(arid1a_sig, arid1b_sig, arid2_sig) 
all_summ <- all_sig %>% group_by(arid, direct) %>% summarise(count = n() ) %>% 
  mutate(percent = count/sum(count) )

all_summ_direct <- all_summ %>% filter(direct == TRUE)
g <- ggplot(all_summ_direct, aes(x = arid, y = percent*100) )
g <- g + geom_bar(stat = 'identity', color = 'grey30', show_guide = FALSE, fill = 'grey50') 
g <- g + theme_figure() + xlab('') + ylab('Percent of differentially expressed\ngenes that are nearest gene to Arid bound region') 
g <- g + coord_fixed(1/10) + theme(axis.text.x = element_text(size = 24, color = 'grey20'), axis.ticks.x = element_blank() )
g
ggsave('output/plots/number_of_diffexpgenes_directtargets.pdf', g)

all_sig %>% 
  filter(direct == TRUE) %>% 
  select(-direct) -> direct_frame

write.table(direct_frame, 'output/diffexp/tables/direct_targets.csv', sep = ',', col.names =T, row.names =F)


