#get per gene exon lengths - for use in normalization of rpkm 
library(data.table) 
library(dplyr) 

extract_gene_name <- function(string) { 
  #string found in Field 9 of the gencode GTF files
  gene_name <- unlist(str_split(string, ';') )
  x <- grep('gene_name', gene_name, perl=T, value=T) %>% 
    str_replace('gene_name', '') %>% 
    str_replace_all('\"', '') %>% 
    str_trim(side = 'both') 
  return(x) 
}

gencode <- fread('~/annotations/gencode.v16.annotation.gtf')
gencode$gene_name <- extract_gene_name(gencode$V9) 
output <- gencode %>% group_by(gene_name) %>% filter(filter = V3 == 'exon') %>% 
                  mutate(length=(V5-V4) ) %>% summarise(exon_length=sum(length) ) 
write.table(output, '~/annotations/gencode.v16.gene_lengths') 
