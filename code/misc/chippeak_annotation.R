#test out chippeakanno library for peak annotation 
library(GenomicRanges) 
library(stringr) 
# below code requires these packages

# functions
# -----------------------------------------------------------------------------
bed2gr <- function(filepath, header=F, sep='\t', stranded=T) { 
  if (!file.exists(filepath) ) { 
    warning('File path does not exist') 
    return(NULL) 
  }
  dat <- read.table(filepath, header=header, sep=sep, stringsAsFactors=FALSE) 
  if (ncol(dat) < 4) {
    dat$V4 <- rep('NONE')
    dat$V5 <- rep(0)
    dat$V6 <- rep('*')
  }
  if (stranded == FALSE) { 
    dat$V6 <- rep('*') 
  } 

  gr <- with(dat, GRanges(seqnames=V1, IRanges(start=V2, end=V3), strand=V6, name=V4) )
  return(gr) 
}

df2gr <- function(df) { 
  # strips df to basic columns (chrom, start, end) and returns a grange object
  require(dplyr)
  colnames(df)[1:3] <- c('a', 'b', 'c') 
  df %<>% select(a, b, c) %>% GRanges(seqnames = a, IRanges(b,c), strand = '*')
  return(df)

}

# -----------------------------------------------------------------------------
extract_gene_name <- function(string) { 
  #string found in Field 10 of the gencode GTF files
  gene_name <- unlist(str_split(string, ';') )
  x <- grep('gene_name', gene_name, perl=T, value=T) %>% 
        str_replace('gene_name', '') %>% 
        str_replace_all('\"', '') %>% 
        str_trim(side = 'both') 
  return(x) 
}

# -----------------------------------------------------------------------------
gtf2gr <- function(filepath, use.data.table=T) { 
  #if the gtf file is the gencode one it is very large, best to use data.table to read it
  #this function is pretty slow so be patient
  if (use.data.table==T) { 
    library(data.table)
    dat <- fread(filepath) 
  } else { 
    dat <- read.table(filepath, comment.char = "#", sep='\t', header=F) 
  }
    dat$name <- extract_gene_name(dat$V9)
  gr <- with(dat, GRanges(seqnames=V1, IRanges(start=V4, end=V5), strand=V7, type=V3, name) )
  dat <- NULL
  gr <- gr[!gr$type %in% c('Selnocysteine', 'start_codon', 'stop_codon')]
  return(gr) 
}

# -----------------------------------------------------------------------------
annotatePeaks <- function(peaks, annotation, window, peak_point = 'center')  { 
  # takes: 
  #  GRanges object of chipseq peaks ()
  #  GRanges object of annotation data (ie gencode genes) 
  #  window controls the range around the tss (default 2kb) 
  # peak_point controls which portion of the peak is checked for overlaps (default 'center') 
  #       also accepts 'all' but this may not behave as expected and is not tested
  # output is the peaks_table with added columns for featureOverlap,
  # -------------------------------------------------------------------------------
  # define new peak_gr based on peak_point
  if (peak_point == 'center') { 
    peaks <- with(peaks, GRanges(seqnames=seqnames,
                                  IRanges(start=floor(start+end)/2, width =2), 
                                  strand = strand, name = name ) )  
  }
  #hierarchy of overlaps <- tss, exon, gene, distal
  #check if peak is in window at TSS
  tss_overlaps <- findNearestFeatureType(peaks, annotation, 'tss', window)
  if (length(tss_overlaps) > 0 ) { 
    tss_overlaps$type <- rep('Promoter') 
    #at each step need to remove those that have a feature annotation 
    hits <- overlapsAny(peaks, tss_overlaps) 
    reduced <- peaks[!hits]
  } else { reduced <- peaks } #just to allow rest of code to work from peak `reduced` object. 
  #check if peak is in exon but not in promoter
  exon_overlaps <- findNearestFeatureType(reduced, annotation, 'exon')
  if (length(exon_overlaps) > 0 ) { 
    exon_overlaps$type <- rep('Exon')
    hits <- overlapsAny(reduced, exon_overlaps)
    reduced <- reduced[!hits]
  }
  #check if peak is in gene but not in exon or promoter - i.e. intron
  intron_overlaps <- findNearestFeatureType(reduced, annotation, 'transcript' )
  if (length(intron_overlaps) > 0 ) {
    intron_overlaps$type <- rep('Intron') 
    hits <- overlapsAny(reduced, intron_overlaps) 
    reduced <- reduced[!hits]
  }
    # if in none of those places call it 'distal'
  distal_overlaps <- reduced
  distal_overlaps$type <- rep('Distal') 
  # combine all the subdata frames and return
  combined <- c(tss_overlaps, exon_overlaps, intron_overlaps, distal_overlaps) 
  #this will break if one of these has no rows
  return(combined) 
}

# -----------------------------------------------------------------------------
getTSS <- function(granges, strand=T, window =1000) { 
  #given a GRange object, return the window aroudn the site of the start site, 
  # true genes should be stranded, 
  #if strand=F will return the 1bp site around the left most coordinate  
  if (strand == T) {
    tss_pos <- granges[strand(granges) == '+',] 
    tss_neg <- granges[strand(granges) == '-', ] 
    newrange_pos <- with(tss_pos, 
                         GRanges(seqnames = seqnames, 
                                 IRanges(start = start - window, 
                                         end = start + window),
                                 strand = strand, name = name) )
    
    newrange_neg <- with(tss_neg ,
                         GRanges(seqnames = seqnames,
                                 IRanges(start = end - window,
                                         end = end + window ), 
                                 strand = strand, name = name) )
    
    adjusted <- c(newrange_pos, newrange_neg)
  } 
  else { 
    adjusted <- with(granges, GRanges(seqname = seqnames, 
                                      IRanges(start = start-window,
                                              end = start + window), 
                                      strand = strand, 
                                      name = name) ) 
                     
  }
  return(adjusted) 
}

# -----------------------------------------------------------------------------
findNearestFeatureType <- function(peaks, annotation, feature, window = 2000) { 
  # takes: 
  # GRanges object of chipseq peaks
  # GRanges object of annotation data (ie gencode genes) 
  # requires a feature to be named; one of (tss | exon | )
  # Returns the chipseq peaks table with added columns for nearest
  if (feature == 'tss') { 
    tss <- annotation[annotation$type == 'transcript', ]  
    adjusted <- getTSS(tss, window = window)
  } 
  else { 
    adjusted <- annotation[annotation$type == feature,]
    #positive vs negative does not matter
  }
  hits <- overlapsAny(peaks, adjusted) 
  return(peaks[hits]) 
}

# -----------------------------------------------------------------------------
distanceToNearestTSS <- function(peaks, annotation, get = TRUE, featurepoint='center') { 
  # takes: 
  # Granges object of peaks
  #  valid metadata is only name - can use that to match up with other function outputs
  #                               
  # Granges object of annotation data (ie gencode genes) 
  #     - one column should be type to select tss if get = T
  #     - if get = F will use whole annotation file
  # featurepoint is whether to measure from center of peak or from closest end.
  
  if (get == TRUE) {
    annotation <- annotation[annotation$type == 'transcript', ]
    tss <- getTSS(annotation, strand = T, window = 1) # only get the point of TSS
  } else {
    tss <- annotation 
  } 
  
  peaks_adjusted <- with(peaks, GRanges(seqnames = seqnames, 
                                        IRanges(start = (start + end)/2, width =2), 
                                        name = name ) ) 
  
  d <- distanceToNearest(peaks_adjusted, tss)
  peaks$distance <- elementMetadata(d)$distance
  nearestGene <- tss[subjectHits(d), ] 
  peaks$nearestTSS <- nearestGene$name
  return(peaks)   
}

# -----------------------------------------------------------------------------


