# Function to return all genes within a given window of a peak or feature
# -----------------------------------------------------------------------
# 
library(GenomicRanges)

genes_in_window <- function(peaks, genes, window = 5e4, promoter =TRUE){ 
  # both should be GRange objects # and genes must have a name field
  # any object in genes that overlaps any part of the area within +/- window/2 will be counted
  # return is a list of unique genes within the window size of peaks
  peak_win <- with(peaks, GRanges(seqnames, IRanges(
                                            start = ((start+end)/2)-(window/2), 
                                            end = ((start+end+1)/2)+(window/2)), 
                                  strand = '*') 
  ) 
  if (promoter == TRUE) { 
    tss_pos <- genes[strand(genes) == '+',] 
    tss_neg <- genes[strand(genes) == '-', ] 
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
    adjusted <- genes
  }
  maps <- findOverlaps(peak_win, genes, select = 'all')
  p <- peaks[queryHits(maps), ]
  g <- genes[subjectHits(maps), ]
  df <- data.frame( chr = seqnames(p), start = start(p), end = end(p), name = p$name, gene_ids  = g$name) 
  return(df)
}

return_expression <- function(genes, expression_df) { 
  # given a list of genes and an expression matrix return mean of expression values
  # expression_df should have a rowname column (made by df %>% add_rownames() )
  df <- expression_df[expression_df$rowname %in% unlist(genes), ] 
  return(mean(df$log2FoldChange, na.rm = T) )
}
