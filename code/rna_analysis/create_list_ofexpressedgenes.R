library(dplyr)
library(ggplot2)
# this script is a helper to make a list of genes expressed vs not that is suitable for use 
# with deep tools
# concatenating the two files will allow one to look at coverage at expressed vs not expressed genes
counts <- read.csv('output/diffexp/tables/count_data_rpkm.csv')
summarised_counts <- counts %>% group_by(gene) %>% summarise(mean = mean(Value, na.rm =T) ) %>%
  arrange(desc(mean) ) %>% 
  filter(!grepl('^MT-.*', gene))
summarised_counts %>% ggplot(aes(x=log2(mean))) + geom_histogram(binwidth = 1) 

all_genes <- read.table('~/annotations/gencode.v16.annotation.bed', header=F, sep = '\t')
all_genes_orderd <- all_genes %>% full_join(summarised_counts, by = c('V4'='gene') ) %>% arrange(desc(mean) )



arid1a_sig <- read.table('~/genes_arid1a_mat', sep ='\t', skip = 1)

expressed_genes <- all_genes_orderd %>% 
  filter(!is.na(mean) )
notexp_genes <- all_genes_orderd %>% 
  filter(is.na(mean) )
write.table(expressed_genes, 'output/diffexp/tables/expressed_genes_ordered_by_exp.bed', sep = '\t', row.names =F, col.names =F, quote =F ) 
write('#Expressed Genes', 'output/diffexp/tables/expressed_genes_ordered_by_exp.bed', append =T) 
write.table(notexp_genes, 'output/diffexp/tables/not_expressed_genes.bed', sep = '\t', row.names = F, col.names =F, quote = F) 
write('#Not Expressed', file = 'output/diffexp/tables/not_expressed_genes.bed', append = T)

