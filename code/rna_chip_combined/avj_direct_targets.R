# One reviewer has asked that the cooperative - competition in figure 1 be done only on direct targets. 
# I'd like to leave the original alone (although I'll change the nomenclature to concordant/discordant to not imply direct action)
# Following the ChIP experiments I will include a new analysis only on direct targets for each gene set (these are few) 
# Additionally I'll include an analysis that is not dependent on peak calls (which are flaky) 

# Concordance of Arid expression changes - DIRECT TARGETS ONLY
# ----------------------------------------------------------------------------
datadir <- 'output/diffexp/tables/'
outputdir <- 'figures/'
source('code/util/theme_paper.R')
source('code/function_defs/rnaseq_fig_fun.R')
library(gridExtra)
library(dplyr)
library(ggplot2)
library(tidyr)

arid1a_sig <- read.csv(paste0(datadir, 'arid1a_sig_results.csv') )
arid1b_sig <- read.csv(paste0(datadir, 'arid1b_sig_results.csv') )
arid2_sig <- read.csv(paste0(datadir, 'arid2_sig_results.csv') )

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

# Now include a filter for genes that are on the peak lists 
arid1a_peaks <- read.csv('output/macs_peaks/cleaned/arid1a_annotated.csv') 
arid1b_peaks <- read.csv('output/macs_peaks/cleaned/arid1b_annotated.csv') 
arid2_peaks  <- read.csv('output/macs_peaks/cleaned/arid2_annotated.csv')
                         
# overwrite above with only those that are directly bound by called peaks - filter extent is after each line
arid1a_arid1b_genes <- arid1a_arid1b_genes[arid1a_arid1b_genes %in% intersect(arid1a_peaks$nearestTSS, arid1b_peaks$nearestTSS)] # 52 -> 15
arid1a_arid2_genes <- arid1a_arid2_genes[arid1a_arid2_genes %in% intersect(arid1a_peaks$nearestTSS, arid2_peaks$nearestTSS)] # 95 -> 37
arid1b_arid2_genes <- arid1b_arid2_genes[arid1b_arid2_genes %in% intersect(arid1b_peaks$nearestTSS, arid2_peaks$nearestTSS)] # 683 -> 161


pal = rev(c('#fc8d59', '#91bfdb'))

#filter log2fc
arid1a_arid1b_df <- getLFC(arid1a_impt, arid1b_impt, c('ARID1A', 'ARID1B'), arid1a_arid1b_genes )
arid1a_arid1b_df %>% 
  do(tidy(cor.test(~ARID1A + ARID1B, data = .) ) ) %>%
  mutate(rsquare = estimate^2)
a <- ggplot(arid1a_arid1b_df, aes(x = ARID1A, y = ARID1B)) 
a <- a + geom_point(aes(fill = ifelse(sign(ARID1A) == sign(ARID1B), 'Concordant', 'Discordant') ), size = 4, shape = 21, color = 'grey30' )
a <- a + theme_figure() + geom_hline(x = 0) + geom_vline(y = 0 )
a <- a + scale_fill_manual(values = pal) 
tmp <- a
a <- a +theme(legend.position = 'none')
a <- a + xlab('ARID1A Log2 Fold Change') + ylab('ARID1B Log2 Fold Change')
a

arid1a_arid2_df <- getLFC(arid1a_impt, arid2_impt, c('ARID1A', 'ARID2'), arid1a_arid2_genes)
arid1a_arid2_df %>% 
  do(tidy(cor.test(~ARID1A + ARID2, data = .) ) ) %>%
  mutate(rsquare = estimate^2)
b <- ggplot(arid1a_arid2_df, aes(x = ARID1A, y = ARID2) )
b <- b + geom_point(aes(fill = ifelse(sign(ARID1A) == sign(ARID2), 'Concordant', 'Discordant') ), size = 4, shape = 21, color = 'grey30')
b <- b + theme_figure() + geom_hline(x =0) + geom_vline(y=0)
b <- b + scale_fill_manual(values = pal) + theme(legend.position = 'none')  
b <- b + xlab('ARID1A Log2 Fold Change') + ylab('ARID2 Log2 Fold Change') 
b

arid1b_arid2_df <- getLFC(arid1b_impt, arid2_impt, c('ARID1B', 'ARID2'), arid1b_arid2_genes) 
c <- ggplot(arid1b_arid2_df, aes(x =  ARID1B, y = ARID2) ) 
c <- c + geom_point(aes(fill = ifelse(sign(ARID1B) == sign(ARID2), 'Concordant', 'Discordant')), size = 4, shape = 21, color = 'grey30') 
c <- c + theme_figure() + geom_hline(x = 0) + geom_vline(y=0)
c <- c + scale_fill_manual(values = pal) + theme(legend.position='none')
c <- c + xlab('ARID1B Log2 Fold Change') + ylab('ARID2 Log2 Fold Change')
c

tmp <- ggplot_gtable(ggplot_build(tmp)) 
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend <- tmp$grobs[[leg]]
pdf('figures/diffexpr_concordance_direct_targets_only.pdf', height = 5, width =15 )
grid.arrange(arrangeGrob(
  arrangeGrob(rectGrob(NA), legend, rectGrob(NA), nrow =1),
  arrangeGrob(a, b, c, nrow = 1), 
  nrow =2, heights = c(1,9) ) ) 
dev.off()


# Get data for all three ARIDs at the ARID1B ARID2 genes 
arid1a_full <- read.csv('output/diffexp/tables/arid1a_full_results.csv')
arid1b_full <- read.csv('output/diffexp/tables/arid1b_full_results.csv') 
arid2_full <- read.csv('output/diffexp/tables/arid2_full_results.csv') 

arid1a_full <- arid1a_full %>% add_rownames() %>% filter(rowname %in% arid1b_arid2_genes) %>% select(rowname, log2FoldChange) 
arid1b_full <- arid1b_full %>%  add_rownames() %>% filter(rowname %in% arid1b_arid2_genes) %>% select(rowname, log2FoldChange) 
arid2_full <- arid2_full %>%  add_rownames() %>% filter(rowname %in% arid1b_arid2_genes) %>% select(rowname, log2FoldChange) 

combined <- arid1a_full %>%
  left_join(arid1b_full, by = 'rowname') %>% 
  left_join(arid2_full, by = 'rowname') 
colnames(combined) <- c('gene', 'ARID1A', 'ARID1B', 'ARID2') 
combined %>%
  arrange(desc(ARID2, ARID1A)) %>%
  mutate(gene = factor(gene, levels = unique(as.character(gene)))) %>%
  gather(arid, value, -gene) %>% 
  ggplot(aes(x = arid, y = gene, fill = value)) + 
  geom_tile() + theme_minimal() + 
  scale_fill_gradient2(low = 'blue2', high = 'red2') + 
  coord_fixed(1/5)

combined %>%
  filter(ARID2 > 0)  
filter(gene %in% rownames(arid1a_sig))
         