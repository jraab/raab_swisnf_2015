# Script to create plots for figure 3 and do statistical testing
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
library(gplots) 
library(readr) 

# Load data sets
# ----------------------------------------------------------------------------
arid1a_all <- read.csv(paste0(datadir, 'arid1a_full_results.csv') )
arid1b_all <- read.csv(paste0(datadir, 'arid1b_full_results.csv') )
arid2_all  <- read.csv(paste0(datadir, 'arid2_full_results.csv') )

arid1a_sig <- read.csv(paste0(datadir, 'arid1a_sig_results.csv') )
arid1b_sig <- read.csv(paste0(datadir, 'arid1b_sig_results.csv') )
arid2_sig <- read.csv(paste0(datadir, 'arid2_sig_results.csv') )

# Function Defs
# -----------------------------------------------------------------------------
test_sig <- function(df) { 
  #just prepass i the filtered data frame
  x <- df %>% 
    select(-arid, -percent) %>% 
    spread(direction, total) %>%
    ungroup() %>% 
    select(Activated, Repressed) 
  print(x)
  return(tidy(chisq.test(x) ) )
}

# Alone vs Jointly by direction
# ----------------------------------------------------------------------------

# Set up data frame with info 
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

write_csv(arid1a_impt, 'output/diffexp/tables/arid1a_diff_genes.csv')
write_csv(arid1b_impt, 'output/diffexp/tables/arid1b_diff_genes.csv') 
write_csv(arid2_impt, 'output/diffexp/tables/arid2_diff_genes.csv') 

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

# Use this to plot alone vs jointly by direction
combined_direction <- combined %>%
  group_by(arid, alone_or_joint, direction) %>%
  summarise(total = n() ) %>% 
  group_by(arid, alone_or_joint) %>% 
  mutate(percent = total/sum(total) * 100  )
pal = c('#fc8d59', '#91bfdb')
h <- ggplot(combined_direction, aes(x= alone_or_joint, y = percent, fill = direction  ) )
h <- h + geom_bar(stat = 'identity') + facet_wrap(~arid, nrow = 1) 
h <- h + geom_bar(stat = 'identity', color = 'grey30', show_guide = F) 
h <- h + theme_minimal() + xlab('') + ylab('Percent of Total Genes') 
h <- h + theme_figure() + theme(strip.text = element_text(size = 20) ) 
h <- h + scale_fill_manual(values = pal) +scale_y_continuous(expand = c(0,0) )
h <- h + theme(legend.direction = 'horizontal', legend.position = 'top', legend.title = element_blank(), 
               panel.grid.major = element_line(color = 'grey10', linetype = 'dashed') ) 
h <- h + coord_fixed(1/20) 
h
ggsave(h, filename = 'figures/fig2A_direction_regulated.pdf', height = 5, width = 5)  
ggsave(h, filename = 'figures/fig2A_direction_regulated.png', dpi = 300) 

#significance test 
combined_direction %>% 
  group_by(arid) %>% 
  do(test_sig(.) )
# Arid1b and Arid2 are very significant (p.values 1.9 e-7, 1.64 e -37 respectively)
# by chisq.test


# Concordance of Arid expression changes 
# ----------------------------------------------------------------------------
arid1a_arid1b_genes <- intersect(arid1a_impt$gene, arid1b_impt$gene) 
arid1a_arid2_genes <- intersect(arid1a_impt$gene, arid2_impt$gene) 
arid1b_arid2_genes <- intersect(arid1b_impt$gene, arid2_impt$gene) 
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
pdf('figures/fig2B_diffexpr_concordance.pdf', height = 5, width =15 )
grid.arrange(arrangeGrob(
  arrangeGrob(rectGrob(NA), legend, rectGrob(NA), nrow =1),
  arrangeGrob(a, b, c, nrow = 1), 
  nrow =2, heights = c(1,9) ) ) 
dev.off()

# Answer to reviewer question - how many ARID1B/ARID2 repressed genes are upregulated after ARID1A
arid1b_up <- row.names(arid1b_sig[arid1b_sig$log2FoldChange > 0,])
arid2_up <- row.names(arid2_sig[arid2_sig$log2FoldChange > 0, ])
arid1a_down <- row.names(arid1a_sig[arid1a_sig$log2FoldChange < 0,] )
arid1a_up <- row.names(arid1a_sig[arid1a_sig$log2FoldChange > 0, ]) 
bothup <- intersect(arid1b_up, arid2_up) 
oppo <- intersect(bothup, arid1a_down)
length(intersect(bothup, arid1a_up) )
arid1a_filt <- arid1a_all[row.names(arid1a_all) %in% bothup, ]
plot(arid1a_filt$log2FoldChange, -log10(arid1a_filt$padj))
points(arid1a_filt[arid1a_filt$padj < 0.1,]$log2FoldChange, -log10(arid1a_filt[arid1a_filt$padj < 0.1,]$padj), col = 'red')

a <- arid1a_all[row.names(arid1a_all) %in% bothup, ]
b <- arid1b_all[row.names(arid1b_all) %in% bothup, ]
c <- arid2_all[row.names(arid2_all) %in% bothup, ]
abc <- data.frame(arid1a=a$log2FoldChange, arid1b=b$log2FoldChange, arid2=c$log2FoldChange) 
pal <- brewer.pal(n = 9, name = 'RdBu') 
heatmap.2(as.matrix(abc), col = pal) 

# add column to alone vs joint data 

