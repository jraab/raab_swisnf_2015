# One reviewer would like a grand total list of all genes changed 
# It should include whether an ARID peak is present (1/0) and the log2FoldChange associated with that gene
# I can also attempt to answer the criticism about figure 1C where I look at alone vs joint genes 
# These genes are not all directly bound - but filtering on binary peak precense underestimates ARID-mediated
# regulation because peak calling is tricky - and peak assignment to genes is trickier. 

library(dplyr) 
library(tidyr) 
library(readr) 
library(gridExtra) 
library(ggplot2) 
source('code/misc/chippeak_annotation.R') # this contains some helper functions 
source('code/util/theme_paper.R')
source('code/function_defs/rnaseq_fig_fun.R')
# Genes
arid1a_genes <- read.csv('output/diffexp/tables/arid1a_sig_results.csv') 
arid1b_genes <- read.csv('output/diffexp/tables/arid1b_sig_results.csv') 
arid2_genes  <- read.csv('output/diffexp/tables/arid2_sig_results.csv') 

# Peaks 
arid1a_peaks <- read.csv('output/macs_peaks/cleaned/arid1a_annotated.csv') 
arid1b_peaks <- read.csv('output/macs_peaks/cleaned/arid1b_annotated.csv') 
arid2_peaks  <- read.csv('output/macs_peaks/cleaned/arid2_annotated.csv') 

# For each genes add column if it is in the peak list. 
arid1a_genes$nearPeak <- ifelse(row.names(arid1a_genes) %in% arid1a_peaks$nearestTSS,  1, 0) 
arid1b_genes$nearPeak <- ifelse(row.names(arid1b_genes) %in% arid1b_peaks$nearestTSS, 1, 0) 
arid2_genes$nearPeak <- ifelse(row.names(arid2_genes) %in% arid2_peaks$nearestTSS, 1, 0) 

# Short cut to name some things as they are in other files
arid1a_sig <- arid1a_genes %>% add_rownames()
rownames(arid1a_sig) <- arid1a_sig$rowname
arid1b_sig <- arid1b_genes %>% add_rownames()
rownames(arid1b_sig) <- arid1b_sig$rowname 
arid2_sig <- arid2_genes  %>% add_rownames()
rownames(arid2_sig) <- arid2_sig$rowname


# Set up data frame with info 
triple <- intersect(row.names(arid1a_sig), intersect(row.names(arid1b_sig), row.names(arid2_sig) ) )

arid1a_alone <- arid1a_sig[!row.names(arid1a_sig) %in% c(row.names(arid1b_sig), row.names(arid2_sig)), ]
arid1b_alone <- arid1b_sig[!row.names(arid1b_sig) %in% c(row.names(arid1a_sig), row.names(arid2_sig)), ]
arid2_alone  <- arid2_sig[!row.names(arid2_sig) %in% c(row.names(arid1a_sig), row.names(arid1b_sig) ), ]

arid1a_jointly <- arid1a_sig[row.names(arid1a_sig) %in% c(row.names(arid1b_sig), row.names(arid2_sig) ), ]
arid1b_jointly <- arid1b_sig[row.names(arid1b_sig) %in% c(row.names(arid1a_sig), row.names(arid2_sig) ), ]
arid2_jointly <- arid2_sig[row.names(arid2_sig) %in% c(row.names(arid1a_sig), row.names(arid1b_sig) ), ]

arid1a_impt <- tabulate_rna_info(arid1a_alone, arid1a_jointly, triple) %>% do(add_number_regulated(.) )  
arid1b_impt <- tabulate_rna_info(arid1b_alone, arid1b_jointly, triple) %>% do(add_number_regulated(.) )  
arid2_impt  <- tabulate_rna_info(arid2_alone, arid2_jointly, triple) %>% do(add_number_regulated(.) )

arid1a_arid1b_genes <- intersect(arid1a_impt$gene, arid1b_impt$gene) 
arid1a_arid2_genes <- intersect(arid1a_impt$gene, arid2_impt$gene) 
arid1b_arid2_genes <- intersect(arid1b_impt$gene, arid2_impt$gene) 

# To address 

# Get Amount of cognate ARID bound at each gene from deeptools output
read_dt <- function(file) { 
  d <- read.table(file, sep = '\t', skip =1, comment.char = '#') 
  x <- d[,7:ncol(d)] 
  row.names(x) <- d[,4]
  colnames(x) <- 1:ncol(x)
  return(x) 
}

# Need to isolate signal just around promoter of genes 
# Files were calculated at -2000:3000genebody:+1000
# Bins were each 10 bp so that means TSS is at bin 200
# Select bins 100:200 which will capture much of signal around gene TSS. 

filter_dt <- function(dt_df) { 
  return(dt_df[,100:250])
  }

arid1a_dt <- read_dt('output/dt_matfiles/genes_arid1a_mat') %>% filter_dt() 
arid1b_dt <- read_dt('output/dt_matfiles/genes_arid1b_mat') %>% filter_dt() 
arid2_dt  <- read_dt('output/dt_matfiles/genes_arid2_mat') %>%  filter_dt() 

# Define 1 number for each peak which can be added to the above
# Median signal around TSS seems reasonable. 
arid1a_med_sig <- arid1a_dt %>% add_rownames() %>%
  gather(bin, value, -rowname) %>% 
  group_by(rowname) %>%  
  summarise(med_signal= median(value, na.rm =T) )

arid1b_med_sig <- arid1b_dt %>% add_rownames() %>%
  gather(bin, value, -rowname) %>% 
  group_by(rowname) %>%  
  summarise(med_signal= median(value, na.rm =T) )

arid2_med_sig <- arid2_dt %>% add_rownames() %>%
  gather(bin, value, -rowname) %>% 
  group_by(rowname) %>%  
  summarise(med_signal= median(value, na.rm =T) )

# Now combine this signal information inot my above information about genes overall. 

arid1a_genes <- arid1a_genes %>% add_rownames() %>% left_join(arid1a_med_sig, by = 'rowname') 
arid1b_genes <- arid1b_genes %>% add_rownames() %>% left_join(arid1b_med_sig, by = 'rowname') 
arid2_genes <- arid2_genes %>% add_rownames() %>% left_join(arid2_med_sig, by = 'rowname') 

# Make lists filtered on the overlap from the AVJ analysis (arid1a + arid1b) 
arid1a_genes_1a_1b <- arid1a_genes %>% filter(rowname %in% arid1a_arid1b_genes) 
arid1b_genes_1a_1b <- arid1b_genes %>% filter(rowname %in% arid1a_arid1b_genes )

arid1a_genes_1a_2 <- arid1a_genes %>% filter(rowname %in% arid1a_arid2_genes) 
arid2_genes_1a_2  <- arid2_genes %>% filter(rowname %in% arid1a_arid2_genes) 

arid1b_genes_1b_2 <- arid1b_genes %>% filter(rowname %in% arid1b_arid2_genes) 
arid2_genes_1b_2  <- arid2_genes  %>% filter(rowname %in% arid1b_arid2_genes) 


# Can now plot these groups to show that for many signficantly changed genes, despite there not being a peak called - 
# There is appreciable signal for the cognate ARID at the genes TSS. This is not bullet proof evidence (since these genes 
# also are not defined as peaks by MACSS), but they do suggest my peak calls are conservative and absence of a 'peak' should not 
# lead to an overconfidence of lack of direct ARID action at a gene. This can be a supplemental figure alongside the direct targeting info
# to demonstrate that point. 
plot_fc <- function(df, alpha = 0.5) { 
  # I put in an alpha so that one can adjust based on number of points
  g<- ggplot(df, aes(x = log2FoldChange, y = log2(med_signal), color = as.factor(nearPeak)) ) + 
    geom_point(size =2, alpha = 0.6) + 
    theme_paper() + 
    theme(panel.grid = element_blank() ) + 
    theme(legend.title = element_text() ) + 
    scale_color_manual(values = c('grey20', 'red2')) + 
    xlab('Log2 FoldChange') + 
    ylab('Median ChIP-seq Signal') +
    guides(color = guide_legend(title = 'Peak Call') )
  return(g) 
  }

legs <- arid1a_genes_1a_1b %>% plot_fc
# I suppresed the legends below - but they can be added back, it was easier to label them in inkscape for the paper

a <- arid1a_genes_1a_1b %>% plot_fc  + theme(legend.position = 'none') + xlab('') + ylab('')
b <- arid1b_genes_1a_1b %>% plot_fc+ theme(legend.position = 'none') + xlab('') + ylab('')

c <- arid1a_genes_1a_2 %>% plot_fc+theme(legend.position = 'none') + xlab('') + ylab('')
d <- arid2_genes_1a_2 %>% plot_fc+ theme(legend.position = 'none') + xlab('')  +ylab('')

e <- arid1b_genes_1b_2 %>% plot_fc (alpha = 0.1) + theme(legend.position = 'none') + xlab('') + ylab('')
f <- arid2_genes_1b_2 %>% plot_fc (alpha = 0.1) + theme(legend.position = 'none') + xlab('') + ylab('')
# now just need to save these plots to a figure. 

pdf('figures/supp_fold_exp_vs_aridsig.pdf') 
grid.arrange(a,b,c,d,e,f, nrow = 3)
dev.off()
