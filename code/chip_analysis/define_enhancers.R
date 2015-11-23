# Create lists of poised vs active enhancers
# The chromatin state maps (chromHMM) are too fragmented

library(dplyr)
library(GenomicRanges) 
library(ggplot2)
source('code/util/chippeak_annotation.R') # need bed2gr function
p300 <- bed2gr('data/external/wgEncodeSydhTfbsHepg2P300sc582IggrabPk.narrowPeak.gz', header = F, sep = '\t', stranded = F) 
k27ac <- bed2gr('data/external/k27ac_hepg2_peaks.narrowPeak', header =F, sep = '\t', stranded = F) 
k4me1 <- bed2gr('data/external/wgEncodeBroadHistoneHepg2H3k04me1StdPk.broadPeak.gz', header =F, sep = '\t', stranded = F) 

# filter only p300 + h3kme1 

all <- p300[overlapsAny(p300, k4me1), ] 

active <- all[overlapsAny(all, k27ac), ]
inactive <- all[!overlapsAny(all, k27ac), ]

# quick check of sizes of these 
summary(width(active)) 
summary(width(inactive) )
#numbers are reasonable

# Need to filter out things that are near to promoters ( I noted a problem when analyzing potential eRNAs)
# It seemed most of my 'hits' which changes in transcription (or maybe any that were transcribed) were
# Within a few kb of a TSS. 
# Hard filter around TSS +/- 5kb
window = 5000

# Get locations of TSS
genes <- bed2gr('output/diffexp/tables/all_genes_ordered_by_exp.bed')
TSS <- getTSS(genes, window = window) 

# Filter out any 'enhancers' within that overlap TSS's within the above window
active_filt <- active[!overlapsAny(active, TSS), ]
inactive_filt <- inactive[!overlapsAny(inactive, TSS) ] 

# write these to a file of enhancers_hepg2
fout <- 'output/misc/hepg2_enhancers_fordt.bed'
write.table(as.data.frame(active_filt), fout, quote = F, sep = '\t', row.names =F, col.names =F) 
write('#Active Enhancers (P300 + H3K4me1 + H3K27ac)', fout, append =T) 
write.table(as.data.frame(inactive_filt), fout, quote= F, sep = '\t', row.names =F, col.names =F, append = T) 
write('#Inactive Enhancers (P300 + H3K4me1)', fout, append = T) 

# Also create a GFF file based on this that can be used with HTSeq-count
active <- as.data.frame(active_filt) 
colnames(active) = c('chr', 'start', 'end', 'score', 'strand', 'x')
inactive <- as.data.frame(inactive_filt) 
colnames(inactive)  = c('chr', 'start', 'end', 'score', 'strand', 'x')


active <- active %>% mutate(name = paste('"', chr, ':', start,'-', end, '"', sep =''), type = 'enhancer', gene_id = 'gene_id', dot = '.', strand2 = '.', src = 'jesse') 
inactive <- inactive %>% mutate(name = paste('"', chr, ':', start, "-", end, '"', sep = ''), type = 'enhancer', gene_id = 'gene_id', dot = '.', strand2 = '.', src = 'jesse')  

active <- active %>% select(chr,src, type, start, end, score, dot, strand2, gene_id, name) %>% 
  mutate(name = paste(gene_id, name, sep =' ')) %>% select(-gene_id) 
inactive <- inactive %>% select(chr,src, type, start, end, score, dot, strand2, gene_id, name, score) %>% 
  mutate(name = paste(gene_id, name, sep =' ')) %>% select(-gene_id) 

write.table(active, 'output/misc/active_enhancers.gff', col.names = F, row.names = F, sep = '\t', quote = F)
write.table(inactive, 'output/misc/inactive_enhancers.gff', col.names = F, row.names =F, sep = '\t', quote = F) 
