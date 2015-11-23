# Figure 8 - Integrate Chip-seq & RNA-seq
library(ggplot2)
library(dplyr) 
library(magrittr)
library(RColorBrewer) 
library(gridExtra)
library(tidyr)
source('code/util/chippeak_annotation.R')
source('code/util/theme_paper.R')
source('code/function_defs/chip_rna_seq_func.R')
pal1 <- rev(c('#ef3b2c', '#fc9272', '#2171b5', '#6baed6'))
cols <- c('#ef3b2c', '#67a9cf', '#7fbf7b',  '#737373') 


# Load Data - requires the following data to be available
# ------------------------------------------------------------------------------
# Peaks
arid1a_peaks <- read.csv('output/macs_peaks/cleaned/arid1a_annotated.csv')
arid1b_peaks <- read.csv('output/macs_peaks/cleaned/arid1b_annotated.csv')
arid2_peaks  <- read.csv('output/macs_peaks/cleaned/arid2_annotated.csv')
multi_peaks  <- read.csv('output/macs_peaks/cleaned/multi_bound_peaks.csv', header =F)
colnames(multi_peaks) <- c('chr', 'start', 'end', 'name', 'nearestTSS', 'distance', 'state')

# make he grange objects
arid1a_gr <- with(arid1a_peaks, GRanges(seqnames = seqnames, IRanges(start, end), strand = '*', name = name))
arid1b_gr <- with(arid1b_peaks, GRanges(seqnames = seqnames, IRanges(start, end), strand = '*', name = name)) 
arid2_gr  <- with(arid2_peaks, GRanges(seqnames = seqnames, IRanges(start, end), strand = '*', name = name) ) 



# Expression  
arid1a_all <- read.csv('output/diffexp/tables/arid1a_full_results.csv') %>% add_rownames()
arid1a_changed <- read.csv('output/diffexp/tables/arid1a_sig_results.csv') %>% add_rownames()
arid1b_all <- read.csv('output/diffexp/tables/arid1b_full_results.csv') %>% add_rownames()
arid1b_changed <- read.csv('output/diffexp/tables/arid1b_sig_results.csv') %>% add_rownames()
arid2_all  <- read.csv('output/diffexp/tables/arid2_full_results.csv') %>% add_rownames()
arid2_changed <- read.csv('output/diffexp/tables/arid2_sig_results.csv') %>% add_rownames()

# This takes a while to load
gencode <- gtf2gr(filepath = 'data/external/gencode.v16.annotation.gtf')
#filter gencode to just the 'gene' objects so that I don't get multiple hits for the same gene/transcript 
gencode2 <- gencode[gencode$type == 'gene', ]

# ------------------------------------------------------------------------------
#function compare expression of genes within adefined region
# this is a terrible function that violates lots of good practices
# it uses environment variables rather than reloading each time
# takes the distance of the neighbordhood and returns the three plots
compare_expression_in_hood <- function(hood_dist) { 
  # get expression info for each peak
  arid1a_p_expr <- neighborhood_expression(arid1a_gr, gencode2, arid1a_all, hood_dist) 
  arid1b_p_expr <- neighborhood_expression(arid1b_gr, gencode2, arid1b_all, hood_dist) 
  arid2_p_expr  <- neighborhood_expression(arid2_gr, gencode2, arid2_all, hood_dist)
  
  arid1a_p_rando <- pick_random_hoods(arid1a_gr, arid1a_all, arid1a_p_expr)
  arid1b_p_rando <- pick_random_hoods(arid1b_gr, arid1b_all, arid1b_p_expr) 
  arid2_p_rando  <- pick_random_hoods(arid2_gr, arid2_all, arid2_p_expr) 
  
  #merge expression into original frame 
  arid1a_peaks_plus <- arid1a_peaks %>% 
    full_join(arid1a_p_expr, by = 'name') %>% 
    tbl_df() %>% 
    mutate(arid = 'Arid1a', distance = hood_dist, random='Arid' )
  
  arid1b_peaks_plus <- arid1b_peaks %>% 
    full_join(arid1b_p_expr, by = 'name') %>% 
    tbl_df() %>% 
    mutate(arid = 'Arid1b', distance = hood_dist, random = 'Arid') 
  
  arid2_peaks_plus <- arid2_peaks %>% 
    full_join(arid2_p_expr, by = 'name') %>% 
    tbl_df() %>% 
    mutate(arid = 'Arid2', distance = hood_dist, random = 'Arid') 
  
  arid1a_rando_plus <- arid1a_peaks %>% 
    full_join(arid1a_p_rando, by = 'name') %>% 
    tbl_df() %>% 
    mutate(arid = 'Arid1a', distance = hood_dist, random = "Random")
  
  arid1b_rando_plus <- arid1b_peaks %>% 
    full_join(arid1b_p_rando, by = 'name') %>%
    tbl_df() %>% 
    mutate(arid = 'Arid1b', distance = hood_dist, random = "Random") 
  
  arid2_rando_plus  <- arid2_peaks %>% 
    full_join(arid2_p_rando, by = 'name') %>% 
    tbl_df() %>% 
    mutate(arid = 'Arid2', distance = hood_dist, random = 'Random') 
  
  all_arids_forhist <- rbind(arid1a_peaks_plus, arid1b_peaks_plus, arid2_peaks_plus,
                             arid1a_rando_plus, arid1b_rando_plus, arid2_rando_plus) 
  
  return(all_arids_forhist)  
} 

# Figure 8A - Total fraction of genes associated with peaks
# ------------------------------------------------------------------------------
arid1a_sig <- arid1a_changed %>%  mutate(direct = rowname %in% arid1a_peaks$nearestTSS, arid = 'Arid1a') 
arid1b_sig <- arid1b_changed %<>% mutate(direct = rowname %in% arid1b_peaks$nearestTSS, arid = 'Arid1b') 
arid2_sig <- arid2_changed %>%  mutate(direct = rowname %in% arid2_peaks$nearestTSS, arid = 'Arid2') 
all_sig <- rbind(arid1a_sig, arid1b_sig, arid2_sig) 
all_summ <- all_sig %>% group_by(arid, direct) %>% summarise(count = n() ) %>% 
  mutate(percent = count/sum(count) )

all_summ_direct <- all_summ %>% filter(direct == TRUE) %>% mutate(arid2 = toupper(arid) ) 
g <- ggplot(all_summ_direct, aes(x = arid2, y = percent*100) )
g <- g + geom_bar(stat = 'identity', color = 'grey30', show_guide = FALSE, fill = 'grey50') 
g <- g + theme_figure() + xlab('') + ylab('Percent of DE Genes') 
g <- g + theme(axis.text.x = element_text(size = 24, color = 'grey20'), axis.ticks.x = element_blank() )
g
ggsave('figures/Fig_7A_number_of_diffexpgenes_directtargets.pdf', g)



hg_nearest(arid1a_changed$rowname, unique(arid1a_peaks$nearestTSS), arid1a_all$rowname ) 
# arid1a p = 0.03
hg_nearest(arid1b_changed$rowname, unique(arid1b_peaks$nearestTSS), arid1b_all$rowname ) 
# arid1b p = 0.67
hg_nearest(arid2_changed$rowname, unique(arid2_peaks$nearestTSS), arid2_all$rowname )
# arid2 p = 3.4e-16


# Figure S8-1 - Single compared to multiple peaks associated with more 
#                gene expression changes (no) 
# ------------------------------------------------------------------------------
arid1a_nearest <- arid1a_peaks %>% check_nearest(arid1a_all)  %>% mutate(arid = 'ARID1A') 
arid1b_nearest <- arid1b_peaks %>% check_nearest(arid1b_all) %>% mutate(arid = 'ARID1B') 
arid2_nearest <- arid2_peaks %>% check_nearest(arid2_all) %>% mutate(arid = 'ARID2') 
all_nearest <- rbind(arid1a_nearest, arid1b_nearest, arid2_nearest) 
all_near_plot <- all_nearest %>%
  ggplot(aes(x= svm, y = percent, fill = isSig)) + geom_bar(position = 'stack', stat = 'identity') + 
  geom_bar(position = 'stack', stat = 'identity', color = 'grey20', show_guide =F) + 
  scale_fill_manual(values = c('grey30', 'steelblue') ) + 
  theme_paper() + facet_wrap(~arid, nrow = 1) + ylab('Percent of Peaks') + xlab('') 
all_near_plot
ggsave('figures/supplemental/S7_1.pdf')

# Plot the change in expression for varying neighborhood sizes (Figure 8B)
# ------------------------------------------------------------------------------
# make plots for each distance
p1e4 <- compare_expression_in_hood(1e4)
p2e4 <- compare_expression_in_hood(2e4) 
p5e4 <- compare_expression_in_hood(5e4)
p1e5 <- compare_expression_in_hood(1e5) 
p2e5 <- compare_expression_in_hood(2e5) 
p5e5 <- compare_expression_in_hood(5e5)
all_distances <- rbind(p1e4, p2e4, p5e4, p2e5, p5e5) 
all_distances %>% ggplot(aes(x = mean_hood_expr, color = svm, fill = svm)) + 
  geom_density(alpha = 0.4) + 
  scale_fill_manual(values = c('grey30', 'steelblue') ) + 
  scale_color_manual(values = c('grey30', 'steelblue') ) + 
  theme_paper() + theme(panel.grid.major.x = element_blank(), 
                        axis.ticks.y = element_blank(), 
                        panel.grid.major.y = element_line(color = 'grey70') )+ 
  xlab('Log2 Fold Change') + ylab('Density') + facet_grid(arid ~ distance) 

dist_plot <- all_distances %>% 
  mutate(distance = distance/1000, arid = toupper(arid) ) %>%
  ggplot(aes(x = as.factor(distance), y = abs(mean_hood_expr), fill = random, color = random) ) + 
  geom_boxplot(notches = T, width = 0.9, color = 'grey20', position = position_dodge(0.9)) + 
  facet_wrap(~arid, nrow =1) + theme_paper() + scale_fill_manual(values = c('steelblue', 'grey30') ) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme(axis.text.x = element_text(color = 'grey10'), 
        axis.ticks = element_blank(), 
        panel.background = element_rect(fill = NA, color = 'grey30')) +
  xlab('Distance Surrounding Peak (kb)') + ylab('Abs Mean Expression\nof Genes in Neighborhood')
ggsave('figures/Fig_7B_dist_plot.pdf', dist_plot, width = 12, height = 4) 

all_distances %>%
  group_by(arid, distance) %>%
  do(tidy(wilcox.test(abs(mean_hood_expr) ~ random, data = .) ))


# Figure 8 Supplement 2 - Comparison of single vs mutliple bound reginos 
#                         # log2 fold change of genes within neighborhood
# -----------------------------------------------------------------------------
aridhist <- all_distances %>% 
  mutate(arid = toupper(arid) ) %>%
  ggplot(aes(x=mean_hood_expr, color = svm, fill = svm) ) +
  geom_density(alpha = 0.25) + 
  scale_fill_manual(values = c('#ef8a62', '#67a9cf') ) + 
  scale_color_manual(values = c('#ef8a62', '#67a9cf')) + 
  theme_paper() + theme(panel.grid.major.x = element_blank(), 
                        axis.ticks.y = element_blank(), 
                        panel.grid.major.y = element_line(color = 'grey70') ) + 
  xlab('') + ylab('') +facet_grid(arid~distance) 
aridhist
ggsave('figures/supplemental/S7_2_svm_distances.pdf') 

# Figure 8 C
# ------------------------------------------------------------------------------
# Stringent set of peaks that are associated to genes regulated and bound by both
# Arid1b and Arid2 (does not require that both bind same place or overlapping)

# Select genes for analysis
rep1b <- arid1b_changed %>% filter(log2FoldChange > 0) %>% select(rowname, log2FoldChange, padj, baseMean) 
rep2 <- arid2_changed %>% filter(log2FoldChange > 0) %>% select(rowname, log2FoldChange, padj, baseMean) 
repressed_in_both <- inner_join(rep1b, rep2, by = 'rowname', x_suffix ='Arid1b', y_suffix = 'Arid2')  
act1b <- arid1b_changed %>% filter(log2FoldChange < 0) %>% select(rowname, log2FoldChange, padj, baseMean) 
act2  <- arid2_changed %>% filter(log2FoldChange < 0 ) %>% select(rowname, log2FoldChange, padj, baseMean)
act_in_both <- inner_join(act1b, act2, by = 'rowname', x_suffix = 'Arid1b', y_suffix = 'Arid2' )

arid1b_peaks_repr_direct <- arid1b_peaks %>% 
  filter(nearestTSS %in% repressed_in_both$rowname) %>% 
  mutate(arid = 'ARID1B', direction = 'Repressed') 

arid2_peaks_repr_direct <- arid2_peaks %>% 
  filter(nearestTSS %in% repressed_in_both$rowname) %>%
  mutate(arid ='ARID2', direction = 'Repressed') 

arid1b_peaks_act_direct <- arid1b_peaks %>%
  filter(nearestTSS %in% act_in_both$rowname) %>%
  mutate(arid = 'ARID1B', direction = 'Activated') 

arid2_peaks_act_direct <- arid2_peaks %>%
  filter(nearestTSS %in% act_in_both$rowname) %>%
  mutate(arid = 'ARID2', direction = 'Activated') 
combined <- rbind(arid1b_peaks_repr_direct, 
                  arid2_peaks_repr_direct, 
                  arid1b_peaks_act_direct, 
                  arid2_peaks_act_direct) 

# What fraction of genes directly regulated and repressed by both are actually directly repressed 
direct_repr_both <- intersect(arid1b_peaks_repr_direct$nearestTSS, arid2_peaks_repr_direct$nearestTSS) 
total_repr_both <- union(arid1b_peaks_repr_direct$nearestTSS, arid2_peaks_repr_direct$nearestTSS) 
length(direct_repr_both)/length(total_repr_both) 
# 42% overlap between directly repressed targets (114/271 total are directly targeted by both both) 

#What fraction of genes directly regulated and activated by both are directly activated
direct_act_both <- intersect(arid1b_peaks_act_direct$nearestTSS, arid2_peaks_act_direct$nearestTSS) 
total_act_both <- union(arid1b_peaks_act_direct$nearestTSS, arid2_peaks_act_direct$nearestTSS) 
length(direct_act_both)/length(total_act_both) 
# 51% of the directly activated targets are shared between both (44/86 genes) 

summary_combined <- combined %>% 
  filter(nearestTSS %in% c(direct_act_both, direct_repr_both) ) %>% 
  group_by(arid, direction, type) %>% 
  summarise(count = n() ) %>% 
  mutate(percent = count/sum(count)) %>%
  mutate(type = factor(type, levels = c('Promoter', 'Exon', 'Distal', 'Intron')))
a1b_v_2_plt <- summary_combined %>% 
  ggplot(aes(x = direction, y = percent, fill = type, order = type)) + 
  geom_bar(position = 'stack', stat = 'identity') + 
  geom_bar(position = 'stack', stat = 'identity', color = 'grey20', show_guide =F) + 
  facet_wrap(~arid, nrow = 1)  + theme_paper() + 
  theme(axis.text = element_text(color = 'grey20')) + xlab('') + ylab('Percent') + 
  scale_fill_manual(values = pal1)
a1b_v_2_plt
ggsave('figures/Fig_7C_1bV2_directreg.pdf', a1b_v_2_plt) 

# Signfiicantly different
summary_combined %>% 
  filter(arid == 'ARID1B') %>%
  select(-percent) %>% 
  spread(direction, count)  %>%  
  select(-type, -arid) %>%  
  do(tidy(chisq.test(cbind(.$Activated, .$Repressed)) ) ) 

# P = 0.08 not diff 
summary_combined %>% 
  filter(arid == 'ARID2') %>% 
  select(-percent) %>%
  spread(direction, count) %>% print() %>% 
  select(-type, -arid) %>% 
  do(tidy(chisq.test(cbind(.$Activated, .$Repressed)))) 


# Supplemental Figure 8D
# -----------------------------------------------------------------------------
expr <- read.csv('output/diffexp/tables/count_data_rpkm.csv') %>% 
  filter(Sample == 'NTG') %>% 
  group_by(gene) %>%
  summarise(Value = mean(Value, na.rm =T) )

arid1a_expr <- arid1a_peaks %>% tbl_df() %>%
  left_join(expr, by =c('nearestTSS' = 'gene') )%>%
  mutate(Value =ifelse(is.na(Value), 0, Value), 
         Group = 'ARID1A') 

arid1b_expr <- arid1b_peaks %>% tbl_df() %>%
  left_join(expr, by = c('nearestTSS' = 'gene') ) %>%
  mutate(Value = ifelse(is.na(Value), 0, Value), 
         Group = 'ARID1B') 

arid2_expr <- arid2_peaks %>% tbl_df() %>%
  left_join(expr, by = c('nearestTSS' ='gene') ) %>%
  mutate(Value = ifelse(is.na(Value), 0, Value), 
         Group = 'ARID2') 


tot_expr <- rbind(arid1a_expr, arid1b_expr, arid2_expr) 
tot_expr_plt <- tot_expr %>% 
  ggplot(aes(x=type, y = log2(Value), fill = svm)) + 
  geom_boxplot(notch=TRUE, position = position_dodge(0.9)) + 
  facet_wrap(~Group, nrow  = 1) +theme_paper() + 
  scale_fill_manual(values = c('', 'steelblue') ) + xlab('') + 
  ylab('Log2 RPKM') +
  theme(panel.background = element_rect(color = 'grey30', fill=NA), 
        axis.text.x = element_text(color = 'grey30', angle = 30, hjust = 1, vjust=1.2), 
        axis.ticks = element_blank() )

tot_expr_plt
ggsave('figures/supplemental/S7_3_total_expr.pdf', tot_expr_plt) 


# calculate p values by wilcox.test for each group pair
tot_expr %>%
  group_by(Group, type) %>%
  summarise(x = as.numeric(pairwise.wilcox.test(log2(Value), svm)$p.value) )
