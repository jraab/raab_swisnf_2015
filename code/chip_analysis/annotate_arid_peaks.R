# Process ChIP-seq Peaks - call from pipeline.chipseq_analysis.sh on codon
# -----------------------------------------------------------------------------
# Take the ChIP-seq peaks called by macs and cleaned up against the blacklist
# Annotate with Genomic Feature
# Assign distance and gene name of nearest transcript 
# Relies on files in helper_scripts/code/annotate_chip_peaks.R
# 
# -----------------------------------------------------------------------------
library(dplyr)
library(tidyr)
# Hard coded locations - fix on your copy of this
#proj_home <- '/magnuson-lab/jraab/analysis/swi_snf_final/' # codon
chip_peaks <-  'output/macs_peaks/cleaned/'
gencode_annotations <- 'data/external/gencode.v16.annotation.gtf'
output_dir <-  'output/macs_peaks/cleaned/'

# Get helper scripts - can clone this repo to $HOME 
#these both need library(GenomicRanges)
source('code/util/chippeak_annotation.R')
source('code/util/merge_grange.R')

# load arid peaks and convert to GRanges objects
arid1a <- bed2gr(paste0(chip_peaks, 'arid1a_peaks_cf.bed'), stranded = F)
arid1b <- bed2gr(paste0(chip_peaks, 'arid1b_peaks_cf.bed'), stranded = F)
arid2  <- bed2gr(paste0(chip_peaks, 'arid2_peaks_cf.bed' ), stranded = F)
snf5   <- bed2gr(paste0(chip_peaks, 'snf5_peaks_cf.bed'), stranded =F)

# load gencode annotations and convert to GRange Object
gencode <- gtf2gr(gencode_annotations)
chromHMM <- bed2gr('output/misc/broadchromhmm_hepg2.bed', stranded =F) # From ENCODE
# see Ernst 2012
# ----------------------------------------------------------------------------
# annotate peaks with genomic feature
window_size <- 1000 # around TSS 
arid1a_annotated <- annotatePeaks(arid1a, gencode, window = window_size) 
arid1b_annotated <- annotatePeaks(arid1b, gencode, window = window_size) 
arid2_annotated  <- annotatePeaks(arid2,  gencode, window = window_size) 
snf5_annotated   <- annotatePeaks(snf5, gencode, window = window_size) 

# annotate peaks with distance to nearest TSS and gene name
arid1a_distance <- distanceToNearestTSS(arid1a, gencode) 
arid1b_distance <- distanceToNearestTSS(arid1b, gencode) 
arid2_distance  <- distanceToNearestTSS(arid2,  gencode) 
snf5_distance   <- distanceToNearestTSS(snf5, gencode)

# label peaks with single if bound by one arid or multiple if bound by more than one
arid1a_distance$svm <- ifelse(overlapsAny(arid1a_distance, union(arid1b_distance, arid2_distance) ), 'Multiple', 'Single') 
arid1b_distance$svm <- ifelse(overlapsAny(arid1b_distance, union(arid1a_distance, arid2_distance) ), 'Multiple', 'Single')
arid2_distance$svm  <- ifelse(overlapsAny(arid2_distance, union(arid1a_distance, arid1b_distance) ), 'Multiple', 'Single') 
# for snf5 this class is overlaps any arid and does not really make any sense
snf5_distance$svm   <- ifelse(overlapsAny(snf5_distance, union(arid1a_distance, union(arid1b_distance, arid2_distance) ) ), 'Multiple', 'Single' )
                              
#also make a peaks set that is the union of all multiple peaks from the 3 arids
# I do not count snf5 here
multi_output <- c(arid1a_distance[arid1a_distance$svm == 'Multiple', ], 
                  arid1b_distance[arid1b_distance$svm == 'Multiple', ], 
                  arid2_distance[arid2_distance$svm   == 'Multiple', ] )
#reduce to a single range for each 
multi_output <- as.data.frame(reduce(multi_output))
names <- sapply(1:nrow(multi_output), function(x) paste('multi', 'peaks', x, sep ='_') )
multi_output$name <- unlist(names)
multi_output$width <- NULL
ordered_names <- c('seqnames', 'start', 'end', 'name')
multi_output <- multi_output[,ordered_names]

# Assign location and gene information for mutli peaks
# ---------------------------------------------------------------------------
multi_gr <- with(multi_output, GRanges(seqnames = seqnames, 
                                       IRanges(start= start, end = end), 
                                       strand  = '*', name = name) )
multi_annotated <- annotatePeaks(multi_gr, gencode, window = window_size)
multi_distance <- distanceToNearestTSS(multi_gr, gencode) 

# Assign chromHMM class to the peaks 
# ---------------------------------------------------------------------------
arid1a_hmm <- annotateChromState(arid1a, chromHMM)
arid1b_hmm <- annotateChromState(arid1b, chromHMM)
arid2_hmm  <- annotateChromState(arid2, chromHMM)
snf5_hmm   <- annotateChromState(snf5, chromHMM)
multi_hmm  <- annotateChromState(multi_gr, chromHMM) 


#merge HMM and feature annotations

# merge distance and feature annotations 
cols_to_keep <- c('name', 'type', 'nearestTSS', 'distance', 'svm') 
arid1a_output <- merge_grange(arid1a_annotated, arid1a_distance, cols_to_keep = cols_to_keep )  
arid1b_output <- merge_grange(arid1b_annotated, arid1b_distance, cols_to_keep = cols_to_keep ) 
arid2_output  <- merge_grange(arid2_annotated, arid2_distance, cols_to_keep = cols_to_keep ) 
snf5_output  <- merge_grange(snf5_annotated, snf5_distance, cols_to_keep = cols_to_keep)
multi_output <- merge_grange(multi_annotated, multi_distance, cols_to_keep = c('name', 'nearestTSS', 'distance') )

# combine the dataframe above with the original grange and prepare 
arid1a_output  <- merge(as.data.frame(arid1a) , arid1a_output, by = 'name') %>% 
          select(seqnames, start, end, name, type, nearestTSS, distance, svm ) %>% 
          merge(as.data.frame(arid1a_hmm, by = 'name')) %>% 
          select(seqnames, start, end, name, type, nearestTSS, distance, svm, state )
arid1b_output <- merge(as.data.frame(arid1b), arid1b_output, by = 'name') %>% 
          select(seqnames, start, end, name, type, nearestTSS, distance, svm) %>%
          merge(as.data.frame(arid1b_hmm, by = 'name')) %>%
          select(seqnames, start, end, name, type, nearestTSS, distance, svm, state)
arid2_output <- merge(as.data.frame(arid2), arid2_output, by = 'name') %>% 
          select(seqnames, start, end , name, type, nearestTSS, distance, svm) %>%
          merge(as.data.frame(arid2_hmm, by = 'name')) %>%
          select(seqnames, start, end, name, type, nearestTSS, distance, svm, state)

snf5_output <- merge(as.data.frame(snf5), snf5_output, by = 'name') %>%
          select(seqnames, start, end, name, type, nearestTSS, distance, svm) %>%
          merge(as.data.frame(snf5_hmm, by = 'name')) %>%
          select(seqnames, start, end, name, type, nearestTSS, distance, svm, state)

multi_output <- merge(as.data.frame(multi_gr), multi_output, by = 'name') %>% 
          select(seqnames, start, end , name, nearestTSS, distance) %>%
          merge(as.data.frame(multi_hmm, by = 'name')) %>%
          select(seqnames, start, end, name,  nearestTSS, distance, state)

# ----------------------------------------------------------------------------
# save output 
#
write.table(arid1a_output, paste0(output_dir, 'arid1a_annotated.csv'),
            sep=',', col.names = T, row.names = F)
write.table(arid1b_output, paste0(output_dir, 'arid1b_annotated.csv'),
            sep=',', col.names = T, row.names = F) 
write.table(arid2_output, paste0(output_dir, 'arid2_annotated.csv'), 
            sep= ',', col.names =T, row.names =F) 
write.table(snf5_output, paste0(output_dir, 'snf5_annotated.csv'), 
            sep = ',', col.names =T, row.names =F) 
write.table(multi_output, paste0(output_dir, 'multi_bound_peaks.csv'), sep=',', col.names =F, row.names = F ) 
# end

