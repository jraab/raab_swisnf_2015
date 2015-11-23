# Functions for RNAseq figures
# Function Defs
# -----------------------------------------------------------------------------
plotsmear <- function(df, label=NULL, thresh = 0.05, xlabs = FALSE){ 
  # df is full output from DESeq2
  sig <- df %>% filter(padj < thresh)
  df <- df %>% filter(padj > thresh) 
  g <- ggplot(df, aes(x= log(baseMean), y = log2FoldChange) ) 
  g <- g + geom_point() 
  g <- g + geom_point(data = sig, aes(x = log2(baseMean), y = log2FoldChange), color = 'red', size = 1.5 ) 
  if (xlabs == TRUE) { 
    g <- g + xlab('Log2 Expression') 
  }
  else
  {
    g <- g + xlab('') 
  }
  g <- g + ylab(label) 
  g <- g + theme_minimal() + geom_hline(y=0, color = 'grey20', linetype = 'dashed') 
  g <- g + theme(axis.text = element_text(size =14), 
                 axis.title = element_text(size = 16) )
  return(g) 
}

# ---------------------------------------------------------
tabulate_rna_info <- function(alone, jointly, triple){ 
  df2 <- data.frame( gene = c(row.names(alone), row.names(jointly) ), 
                     baseMean = c(alone$baseMean, jointly$baseMean) , 
                     log2FoldChange = c(alone$log2FoldChange, jointly$log2FoldChange), 
                     padj = c(alone$padj, jointly$padj), 
                     alone_or_joint = c(rep('Alone', nrow(alone)), rep('Jointly', nrow(jointly) ) ) )
  df2 <- df2 %>% mutate(direction = ifelse(df2$log2FoldChange < 0, 'Activated', 'Repressed'), 
                        inall = ifelse(df2$gene %in% triple, TRUE, FALSE) ) 
  return(df2) 
}

# --------------------------------------------------------------------------------------
add_number_regulated <- function(df){ 
  df$num_regulated <- ifelse(df$alone_or_joint == 'Jointly' & df$inall == TRUE, 'Three', 'Two') 
  df$num_regulated <- ifelse(df$num_regulated == 'Two' & df$alone_or_joint == 'Alone', 'One', df$num_regulated)
  return(df) 
}

# ----------------------------------------------------------------------------
getLFC <- function(df1, df2, names, selector ) { 
  df1_sel <- df1 %>% filter(gene %in% selector) %>% select(log2FoldChange)
  df2_sel <- df2 %>% filter(gene %in% selector) %>% select(log2FoldChange)
  dfref <- cbind(df1_sel, df2_sel) 
  colnames(dfref) <- names 
  return(dfref) 
}