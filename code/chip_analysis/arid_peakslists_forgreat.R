# Make lists of ARID peaks
library(readr)
library(dplyr)

#
if (!dir.exists('output/macs_peaks/cleaned/for_great/')) { dir.create('output/macs_peaks/cleaned/for_great')}

arid1b <- read_csv('output/macs_peaks/cleaned/arid1b_annotated.csv')
arid1b_single <- arid1b%>% filter(svm == 'Single') %>% select(seqnames, start, end, name)
write.table(arid1b_single, 'output/macs_peaks/cleaned/for_great/arid1b_single.bed', row.names = F, sep = '\t', col.names = F, quote = F)

