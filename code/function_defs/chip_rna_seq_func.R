# Functions for comparing ChIP-seq and RNA-seq

hg_nearest <- function(altered_genes, genes_nearest_peaks, all_genes) { 
  # calculate a p value for whether more altered genes than would be expected by chance are 
  # the nearest gene to a peak of the same arid
  # In each case the genes above should be passed as list
  total <- length(all_genes)
  overlap <- length(intersect(altered_genes, genes_nearest_peaks))
  print(overlap/length(altered_genes) )
  all_altered <- length(altered_genes) 
  all_near <- length(genes_nearest_peaks) 
  p <- phyper(overlap, all_altered, total, all_near, lower.tail = F)
  return(p) 
}

return_expression <- function(genes, expression_df) { 
  # given a list of genes and an expression matrix return mean of expression values
  # expression_df should have a rowname column (made by df %>% add_rownames() )
  df <- expression_df[expression_df$rowname %in% unlist(genes), ] 
  return(mean(df$log2FoldChange, na.rm = T) )
}

check_nearest <- function(peak_df, expr_df) { 
  # Is the nearest gene to a chip-seq peak significantly changed in RNAseq data
  peak_df %>% 
    full_join(expr_df, by = c('nearestTSS' = 'rowname')) %>% 
    tbl_df() %>% 
    select(name, type, nearestTSS, distance, svm, log2FoldChange, pvalue) %>% 
    filter(!is.na(log2FoldChange), !is.na(type) ) %>%
    mutate(isSig = ifelse(pvalue < 0.05, 'Significant', 'Not Significant') ) %>% 
    group_by(svm, isSig) %>% 
    summarise(count = n() ) %>% 
    mutate(percent = 100* count/sum(count) )
}

# function to assign gene expression values for all peaks with all genes in neighborhood. 
neighborhood_expression <- function(peaks_gr, gene_info, full_expression_table, window){ 
  require(dplyr) 
  require(magrittr) #loaded separately for the %<>% operator
  source('code/util/define_genes_in_window.R') # This and other helper functions need packaged
  
  # function takes a GRanges object of peaks and looks up the expression values for all
  # genes within a defined window size
  # expression table should contain all expressed genes for the experiment
  # as implmented returns the mean(log2FoldChange) and full_expression must have that column
  # peaks_gr  = GRanges object of peak locations and names - names will be used in return to facilitate mapping
  # gene_info = GRanges object of gene information - such as the gencode annotations from gtf2gr function
  # full_expression = data.frame containing information about log2FoldExpression for genes in experiment. 
  # window = size of window around peak center to look - overlaps currently defined as any TSS within this window as assigned to peak. 
  df <- genes_in_window(peaks_gr, gene_info, window = window) 
  df %<>% group_by(name) %>% 
    summarise(genes_in_hood = list(as.character(gene_ids))) %>% 
    group_by(name) %>% 
    mutate(mean_hood_expr = return_expression(genes_in_hood, full_expression_table) )
  return(df) 
}

pick_random_hoods <- function(gr, expr, p_expr) {
  # need a way to get random expression for a set of peaks (arid1a, arid1b, arid2) .
  # this function gets called during compare_expression_in_hood and requires the 
  # output of neighborhood_expreesion to give the distirbution of number of gnees
  # I'll randomly select that number of genes from the expr data and return the data
  gene_num_dist <- unlist(lapply(p_expr$genes_in_hood, function(x) length(unlist(x) ) )) 
  num_to_sample <- sample(gene_num_dist, length(gr), replace = T )
  # get a list of genes corresponding to these numbers
  vals <- lapply(num_to_sample, function(x) mean(expr[sample(1:nrow(expr), size = x, replace =F), 'log2FoldChange'] ))
  gr %<>% as.data.frame() %>% tbl_df() %>%
    mutate(mean_hood_expr = unlist(vals), genes_in_hood = rep('NA')) %>%
    select(name, genes_in_hood, mean_hood_expr) 
  return(gr) 
}