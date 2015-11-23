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
merge_grange <- function(gr1, gr2, cols_to_keep) { 
  m <- merge(as.data.frame(gr1), as.data.frame(gr2), by = 'name') 
  m <- m[,cols_to_keep]
  return(m) 
} 

# ---------------------------------------------------------------------------
draw2way <- function(gr1, gr2, cat_names, colors ) { 
  ov <- length(intersect(gr1, gr2) )
  grid.newpage()
  draw.pairwise.venn(length(gr1), length(gr2), cross.area = ov, category = cat_names, 
                     scaled = T, col = 'grey30', fill = colors, alpha = 0.5, 
                     label.col = 'Grey20', cex = 2, lwd = 2, cat.cex = 2, cat.dist = rep(0.025, 2), 
                     cat.fontfamily = rep('Helvetica', 2)) 
}

# -----------------------------------------------------------------------------
draw3way <- function(gr1, gr2, gr3, labels , colors) {
  n12 <- length(intersect(gr1, gr2) )  
  n13 <- length( intersect(gr1, gr3) )
  n23 <- length(intersect(gr2, gr3) )
  n123 <- length(intersect(gr1, intersect(gr2, gr3) ) )
  grid.newpage()
  draw.triple.venn(length(gr1), length(gr2), length(gr3), 
                   n12, n23, n13, n123, 
                   category = labels, col = 'grey30', fill = colors, alpha = 0.5,
                   cex = 2, lwd = 2, cat.cex = 2, cat.dist = rep(0.025, 3), 
                   cat.fontfamily = rep('Helvetica', 3) )
  return()
}

# Functions for metaplotting from hdf5 
# -----------------------------------------------------------------------------
readh5vals <- function(fn) {
  tmp <- h5read(fn, 'coverages', bit64conversion='double')
  vals <- t(tmp$block0_values) 
  row.names(vals) <- tmp$axis1
  return(vals) 
}

readh5vals_fortype <- function(fn) {
  tmp <- h5read(fn, 'coverages', bit64conversion='double')
  vals <- data.frame(t(tmp$block0_values) )
  vals$rowname <- tmp$axis1
  print(dim(vals) )
  return(vals) 
}

filter_arid_coverage <- function(arid_coverage_file, peak_names) { 
  # will return just the coverage for the correct list of peaks
  vals <- readh5vals(arid_coverage_file) 
  return(vals[row.names(vals) %in% peak_names]) 
}

collapsedf <- function(list_dfs) { 
  odf <- data.frame()
  for (i in 1:length(list_dfs) ) { 
    odf <- rbind(odf, list_dfs[[i]])
  }
  return(odf)
}

melt_df_list <- function(df_list) { 
  for (i in 1:length(df_list)){
    df <- df_list[[i]]
    df$Peakname <- row.names(df) 
    mf <- gather(df, 'Position', 'Value', -Peakname) 
  }
  return(output_df) 
}

ci <- function(vals) { 
  return( qnorm(0.975) * sd(vals, na.rm =T) /sqrt(length(vals) ) )
}



plot_four_groups <- function(end, stubs, peaknames) { 
  cols <- c('#ef3b2c', '#67a9cf', '#7fbf7b',  '#737373') 
  vals <- sapply(stubs, function(x) readh5vals(paste0('output/encode_coverages/coverages/', x, end)))
  #vals needs filtered
  
  o <- lapply(seq_along(vals),
              function(x, v) { as.data.frame(v[[x]]) %>% mutate(PeakNames = row.names(.) ) %>% 
                  filter(PeakNames %in% peaknames[[x]]) %>%
                  gather('Position', 'Value', -PeakNames) %>% 
                  mutate(group = as.character(rep(names(v[x])) )  ) %>% 
                  group_by(Position, group) %>% 
                  summarise(avg = mean(Value, na.rm =T), ci = ci(Value) ) %>% 
                  mutate (ymin = avg - ci, ymax = avg +ci )}, v = vals) 
  
  output_df <- collapsedf(o) 
  output_df$Position <- rep(seq(-2500, 2490, by = 10), 4) 
  output_df$group <- sapply(output_df$group, toupper) 
  g <- ggplot(output_df, aes(x = Position, y = avg, color = group, fill = group) )
  g <- g + geom_line(show_guide = T, size = 1.5) 
  g <- g +  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.25, color = NA ) 
  g <- g + guides(color = guide_legend(overide.aes = list(size =4) ) )
  g <- g +theme_minimal() + xlab('') + ylab('IP/Input') 
  g <- g + theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank(), 
                 axis.text = element_text(size = 20), axis.title = element_text(size = 24), 
                 legend.text = element_text(size = 24))  
  g <- g + scale_fill_manual(values = cols) 
  g <- g + scale_color_manual(values = cols) 
  return(g) 
}

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}

plot_arids <- function(end, stubs) { 
  # plot metagene plots for all aird peaks for a given hdf5 files
  cols <- c('#ef3b2c', '#67a9cf', '#7fbf7b',  '#737373') 
  vals <- sapply(stubs, function(x) readh5vals(paste0('output/encode_coverages/coverages/', x, end) ) ) 
  o <- lapply(seq_along(vals), 
              function(x, v) { as.data.frame(v[[x]]) %>%  
                  mutate(PeakNames = row.names(.) ) %>% 
                  gather('Position', 'Value', -PeakNames) %>% 
                  mutate(group = as.character(rep(names(v[x] ) ) ) ) %>% 
                  group_by(Position, group) %>% 
                  summarise(avg = mean(Value, na.rm = T), ci = ci(Value) ) %>% 
                  mutate (ymin = avg - ci, ymax = avg + ci)
              }, v = vals )
  output_df <- collapsedf(o) 
  output_df$Position <- rep(seq(-2500, 2490, by = 10), length(stubs)) # only 3 groups this time
  output_df$group <- sapply(output_df$group, toupper) 
  g <- ggplot(output_df, aes(x = Position/1000, y = avg, color = group, fill = group) )
  g <- g + geom_line(show_guide = T, size = 1.5) 
  g <- g +  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.25, color = NA, show_guide =F ) 
  g <- g + guides(color = guide_legend(overide.aes = list(size =8) ) )
  g <- g +theme_minimal() + xlab('') + ylab('IP/Input') 
  g <- g + theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank(), 
                 axis.text = element_text(size = 20), axis.title = element_text(size = 24), 
                 legend.text = element_text(size = 24))  
  g <- g + scale_fill_manual(values = cols) 
  g <- g + scale_color_manual(values = cols) 
  return(g) 
}


plot_arids_by_type <- function(end, stubs, peak_info) { 
  # same as above but plot the arids distinguised by genomic feature
  
  vals <- lapply(stubs, function(x) readh5vals_fortype(paste0(dirad, x, end) ) ) 
  #vals <- lapply(vals, function(x) as.data.frame(x) )
  #vals <- lapply(vals, function(x) x$rowname <- row.names(x) )
  o <- lapply(seq_along(vals), 
              function(x, v) { v[[x]] %>% tbl_df() %>% 
                  full_join(peak_info[[x]], by = c('rowname'='name')) %>% 
                  select(-nearestTSS, -distance, -svm, -state, -start, -end, -seqnames) %>% 
                  gather('Position', 'Value', -rowname, -type) %>%
                  mutate(group = rep(names(peak_info[x] ) ) )  %>% 
                  filter(!type == 'Exon') %>%
                  group_by(Position, group, type) %>% 
                  summarise(avg = mean(Value, na.rm = T), ci = ci(Value) ) %>% 
                  mutate (ymin = avg - ci, ymax = avg + ci) 
              }, v = vals )
  output_df <- collapsedf(o) 
  
  output_df$Position <- rep(seq(-2500, 2490, by = 10), each = 3, times = 3) # only 3 groups this time
  g <- ggplot(output_df, aes(x = Position, y = avg, color = type , fill = type) )
  g <- g + geom_line(show_guide = F, size = 1.5) +  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.5, color = NA ) 
  g <- g +theme_minimal() + xlab('') + ylab('IP/Input') +theme_paper()
  g <- g + theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank(), 
                 axis.text = element_text(size = 20), axis.title = element_text(size = 24), 
                 legend.text = element_text(size = 24))  
  g <- g + facet_wrap(~group,nrow =1) 
  
  return(g) 
}

getfn <- function(path, rt_patt) { 
  # return list of files that match the pattern
  patt <- paste0('.*_p_',rt_patt,'.*_s_coverage.h5')
  files <- list.files(path = path, patt = patt, full.names = T, ignore.case = T )
  files <- files[!grepl(pattern = 'snf5*', files)] 
  return(files) 
}

getFilenames <- function(arid, endings) { 
  # return list of files that match the pattern
  path <- 'output/encode_coverages/coverages/'
  file_patt <- lapply(endings, function(x) paste0(arid, x) ) 
  files <- unlist(lapply(file_patt, function(x) list.files(path, x, full.names =T, ignore.case =T) ))
  print (files) 
  return(files) 
}
# ------------------------------------------------------------------------------