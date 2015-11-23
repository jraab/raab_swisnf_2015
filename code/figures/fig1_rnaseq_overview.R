# Script needed to create most of Figure 2

datadir <- 'output/diffexp/tables/'
outputdir <- 'figures/'
source('code/util/theme_paper.R')
source('code/function_defs/rnaseq_fig_fun.R')
library(ggplot2) 
library(dplyr) 
library(tidyr)
library(gridExtra)
library(RColorBrewer)
library(broom)

# Load data sets
# ----------------------------------------------------------------------------
arid1a_all <- read.csv(paste0(datadir, 'arid1a_full_results.csv') )
arid1b_all <- read.csv(paste0(datadir, 'arid1b_full_results.csv') )
arid2_all  <- read.csv(paste0(datadir, 'arid2_full_results.csv') )

arid1a_sig <- read.csv(paste0(datadir, 'arid1a_sig_results.csv') )
arid1b_sig <- read.csv(paste0(datadir, 'arid1b_sig_results.csv') )
arid2_sig <- read.csv(paste0(datadir, 'arid2_sig_results.csv') )

# Arid fold change from RNA-seq data
# ------------------------------------------------------------------------------
#check the RNAseq data for each 
groups <- c('ARID1A', 'ARID1B', 'ARID2') 
arid1a_all_rs <- arid1a_all  %>%
  add_rownames() %>% 
  filter(rowname %in% groups) %>% 
  select(rowname, log2FoldChange) %>% 
  mutate(si = 'siARID1A') 

arid1b_all_rs <- arid1b_all %>% 
  add_rownames() %>% 
  filter(rowname %in% groups) %>% 
  select(rowname, log2FoldChange) %>% 
  mutate(si = 'siARID1B')

arid2_all_rs <- arid2_all %>% 
  add_rownames() %>% 
  filter(rowname %in% groups) %>% 
  select(rowname, log2FoldChange) %>% 
  mutate(si = 'siARID2')

combined <- rbind(arid1a_all_rs, arid1b_all_rs, arid2_all_rs) 
extra <- combined[combined$rowname == toupper(gsub('si', '', combined$si) ), ]
cols <-c('#ef3b2c', '#67a9cf', '#7fbf7b') 
g <- combined %>% 
  ggplot(aes(x = rowname, y = log2FoldChange, fill = si) ) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  geom_bar(stat = 'identity', position = 'dodge', color = 'grey30', show_guide = F) + 
  theme_paper() + theme(axis.text = element_text(size = 16, color = 'grey20'), 
                        axis.ticks = element_blank() )+ 
  xlab('Gene Expression') + ylab('Log2 Fold Change') + 
  scale_fill_manual(values = cols) 
g    

ggsave('figures/1b_rnaseq_foldchange_arids.pdf', g ,height = 4, width = 4) 

# Smear plot - function defined above ( Figure 2C )
# ----------------------------------------------------------------------------
arid1a_smear <- plotsmear(arid1a_all, 'ARID1A', 0.1)
arid1b_smear <- plotsmear(arid1b_all, 'ARID1B', 0.05) 
arid2_smear  <- plotsmear(arid2_all, 'ARID2', 0.05, xlabs = TRUE) 
pdf(paste0(outputdir, 'fig1c_smearplot.pdf') ) 
grid.arrange(arid1a_smear, arid1b_smear, arid2_smear, ncol = 1) 
dev.off()

# Number Upregulated/Downregulated 
# Not saved for now.
all_sig <- rbind(arid1a_sig, arid1b_sig, arid2_sig) 
all_sig$dir <- ifelse(sign(all_sig$log2FoldChange) > 0, 'Up Regulated', 'Down Regulated')
all_sig$arid <- c(rep('ARID1A', nrow(arid1a_sig)), 
                  rep('ARID1B', nrow(arid1b_sig)), 
                  rep('ARID2', nrow(arid2_sig)) )
ggplot(all_sig, aes(x=arid, fill = dir) ) + geom_bar(position = 'dodge', stat = 'bin')  + 
  theme_paper()

# Alone vs Jointly - total numbers ( Figure 2D )
# ----------------------------------------------------------------------------
triple <- intersect(row.names(arid1a_sig), intersect(row.names(arid1b_sig), row.names(arid2_sig) ) )

arid1a_alone <- arid1a_sig[!row.names(arid1a_sig) %in% c(row.names(arid1b_sig), row.names(arid2_sig)), ]
arid1b_alone <- arid1b_sig[!row.names(arid1b_sig) %in% c(row.names(arid1a_sig), row.names(arid2_sig)), ]
arid2_alone  <- arid2_sig[!row.names(arid2_sig) %in% c(row.names(arid1a_sig), row.names(arid1b_sig) ), ]

arid1a_jointly <- arid1a_sig[row.names(arid1a_sig) %in% c(row.names(arid1b_sig), row.names(arid2_sig) ), ]
arid1b_jointly <- arid1b_sig[row.names(arid1b_sig) %in% c(row.names(arid1a_sig), row.names(arid2_sig) ), ]
arid2_jointly <- arid2_sig[row.names(arid2_sig) %in% c(row.names(arid1a_sig), row.names(arid1b_sig) ), ]

# build data frames with important info 


arid1a_impt <- tabulate_rna_info(arid1a_alone, arid1a_jointly, triple) %>% do(add_number_regulated(.) )  
arid1b_impt <- tabulate_rna_info(arid1b_alone, arid1b_jointly, triple) %>% do(add_number_regulated(.) )  
arid2_impt  <- tabulate_rna_info(arid2_alone, arid2_jointly, triple) %>% do(add_number_regulated(.) )

arid1a_impt <- arid1a_impt 
arid1b_impt <- arid1b_impt
arid2_impt  <- arid2_impt 

# combined data frames for each arid 
combined <- tbl_df(data.frame(arid = c(rep('ARID1A', nrow(arid1a_impt)), 
                                       rep('ARID1B', nrow(arid1b_impt)), 
                                       rep('ARID2',  nrow(arid2_impt) ) ),
                              num_reg = c(arid1a_impt$num_regulated, arid1b_impt$num_regulated, arid2_impt$num_regulated), 
                              alone_or_joint = c(as.character(arid1a_impt$alone_or_joint), 
                                                 as.character(arid1b_impt$alone_or_joint), 
                                                 as.character(arid2_impt$alone_or_joint) ), 
                              direction = c(arid1a_impt$direction, arid1b_impt$direction, arid2_impt$direction), 
                              log2FC  = c(arid1a_impt$log2FoldChange, arid1b_impt$log2FoldChange, arid2_impt$log2FoldChange) ) )
# plot total numbers 
combined_totals <- combined %>% group_by(arid, num_reg) %>% 
  summarise(total = n() ) %>%
  group_by(arid) %>% 
  mutate(percent = total/sum(total) * 100)

combined_totals$num_reg <- factor(combined_totals$num_reg, levels = c('One', 'Two', 'Three') ) 
combined_totals <- combined_totals %>% arrange(num_reg)
print(combined_totals)

pal <- rev(brewer.pal(3, 'Blues') )
g <- ggplot(combined_totals, aes(x = arid, y = percent, fill = num_reg) )
g <- g + geom_bar(stat = 'identity')
g <- g + geom_bar(stat = 'identity', color = 'grey10', show_guide = F)
g <- g + theme_classic() + xlab('') + ylab('Percent') 
g <- g + theme_paper()
g <- g + scale_fill_manual(values = pal) +scale_y_continuous(expand = c(0,0) )
g
ggsave(filename = paste0(outputdir, 'fig2D_number_regulated.pdf'), g) 

# As a pie chart as requested by a reviewer 
pieversion <- ggplot(combined_totals, aes(x = factor(1), y = percent, fill = num_reg )) +
  geom_bar(width = 1, stat = 'identity', color = 'grey30', show_guide=F) + 
  geom_bar(width = 1, stat = 'identity') + 
  coord_polar(theta = 'y') + 
  facet_wrap(~arid) + 
  scale_fill_manual(values = pal) + 
  theme_minimal()  + xlab('') + ylab('')  + 
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid = element_blank(), 
        legend.title = element_blank(), 
        strip.text = element_text(size = 20) 
        )

ggsave('paper/figures/fig1a_pie.pdf', pieversion)

# Few bits of important information
# ------------------------------------------------------------------------------
# get total number of regulated genes. 
all_unique_genes <- length(unique(c( row.names(arid1a_sig), row.names(arid1b_sig), row.names(arid2_sig)))) 
print(all_unique_genes)

# and perecentage of 'expressed' genes
print(all_unique_genes/length(row.names(arid1a_all)))
# ie. 20%