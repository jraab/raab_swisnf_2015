# Intersection of RNA seq and ChIP-seq data
# -----------------------------------------------------------------------------
library(ggplot2)
library(dplyr) 
library(magrittr)
library(RColorBrewer) 
library(gridExtra)
library(tidyr)
library(broom)
source('/Users/jraab/helper_scripts/code/chippeak_annotation.R')
source('code/util/theme_paper.R')
pal1 <- rev(c('#ef3b2c', '#fc9272', '#2171b5', '#6baed6'))
cols <- c('#ef3b2c', '#67a9cf', '#7fbf7b',  '#737373') 

#source('/Users/jraab/helper_scripts/code/')
# Load Data
# ------------------------------------------------------------------------------
# Peaks
arid1a_peaks <- read.csv('output/macs_peaks/cleaned/arid1a_annotated.csv')
arid1b_peaks <- read.csv('output/macs_peaks/cleaned/arid1b_annotated.csv')
arid2_peaks  <- read.csv('output/macs_peaks/cleaned/arid2_annotated.csv')
multi_peaks  <- read.csv('output/macs_peaks/cleaned/multi_bound_peaks.csv', header =F)
colnames(multi_peaks) <- c('chr', 'start', 'end', 'name', 'nearestTSS', 'distance', 'state')

# Expression  
arid1a_all <- read.csv('output/diffexp/tables/arid1a_full_results.csv') %>% add_rownames()
arid1a_changed <- read.csv('output/diffexp/tables/arid1a_sig_results.csv') %>% add_rownames()
arid1b_all <- read.csv('output/diffexp/tables/arid1b_full_results.csv') %>% add_rownames()
arid1b_changed <- read.csv('output/diffexp/tables/arid1b_sig_results.csv') %>% add_rownames()
arid2_all  <- read.csv('output/diffexp/tables/arid2_full_results.csv') %>% add_rownames()
arid2_changed <- read.csv('output/diffexp/tables/arid2_sig_results.csv') %>% add_rownames()

# Basic Summary 
# -----------------------------------------------------------------------------
arid1a_peaks %>% group_by(type, svm) %>% 
  filter(nearestTSS %in% arid1a_changed$rowname) %>% 
  summarise(count = n())  
 
arid1b_peaks %>% group_by(type, svm) %>% 
  filter(nearestTSS %in% arid1b_changed$rowname) %>% 
  summarise(count = n() )

arid2_peaks %>% group_by(type, svm) %>% 
  filter(nearestTSS %in% arid2_changed$rowname ) %>% 
  summarise(count = n() )

# Look at genes regulated by pairs of Arids
# ---------------------------------------------------------------------------
# Arid1a and Arid2 
g_1a_2 <- arid1a_changed[row.names(arid1a_changed) %in% row.names(arid2_changed),'rowname']  

g_1b_2 <-arid1b_changed[row.names(arid1b_changed) %in% row.names(arid2_changed), 'rowname' ]
# how many direct targets can we id.
arid1b_peaks %>% filter(nearestTSS %in% g_1b_2) %>% select(nearestTSS) %>% unique(.) %>% count() #379 unique genes
arid2_peaks  %>% filter(nearestTSS %in% g_1b_2) %>% select(nearestTSS) %>% unique(.) %>% count() # 583 unique genes
 

# Of 680 genes that are altered in both arid1b and arid2 219(arid1b) and 305(arid2) are the nearest gene to at least 1 arid1b or arid2 peak respectively.
gencode <- gtf2gr(filepath = '/Users/jraab/annotations/gencode.v16.annotation.gtf')
#filter gencode to just the 'gene' objects so that I don't get multiple hits for the same gene/transcript 
gencode2 <- gencode[gencode$type == 'gene', ]

# Are genes that are nearest to Arid peaks more likely to be changed 
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

hg_nearest(arid1a_changed$rowname, unique(arid1a_peaks$nearestTSS), arid1a_all$rowname ) 
# arid1a p = 0.03
hg_nearest(arid1b_changed$rowname, unique(arid1b_peaks$nearestTSS), arid1b_all$rowname ) 
# arid1b p = 0.67
hg_nearest(arid2_changed$rowname, unique(arid2_peaks$nearestTSS), arid2_all$rowname )
# arid2 p = 3.4e-16

mean(abs(arid1a_all[arid1a_all$rowname %in% arid1a_peaks$nearestTSS, 'log2FoldChange']), na.rm =T)
mean(abs(arid1b_all[arid1b_all$rowname %in% arid1b_peaks$nearestTSS, 'log2FoldChange']), na.rm =T) 
mean(abs(arid2_all[arid2_all$rowname %in% arid2_peaks$nearestTSS, 'log2FoldChange']), na.rm =T) 

# Are more genes that are direct targets found for different classes 
# Foreach peak data frame create new data frame 
# with peak name, nearestTSS, svm, type, log2FoldChange, isSig

return_expression <- function(genes, expression_df) { 
  # given a list of genes and an expression matrix return mean of expression values
  # expression_df should have a rowname column (made by df %>% add_rownames() )
  df <- expression_df[expression_df$rowname %in% unlist(genes), ] 
  return(mean(df$log2FoldChange, na.rm = T) )
}
  
check_nearest <- function(peak_df, expr_df) { 
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

arid1a_nearest <- arid1a_peaks %>% check_nearest(arid1a_all)  %>% mutate(arid = 'Arid1a') 
arid1b_nearest <- arid1b_peaks %>% check_nearest(arid1b_all) %>% mutate(arid = 'Arid1b') 
arid2_nearest <- arid2_peaks %>% check_nearest(arid2_all) %>% mutate(arid = 'Arid2') 
all_nearest <- rbind(arid1a_nearest, arid1b_nearest, arid2_nearest) 
all_near_plot <- all_nearest %>%
  ggplot(aes(x= svm, y = percent, fill = isSig)) + geom_bar(position = 'stack', stat = 'identity') + 
    geom_bar(position = 'stack', stat = 'identity', color = 'grey20', show_guide =F) + 
    scale_fill_manual(values = c('grey30', 'steelblue') ) + 
    theme_paper() + facet_wrap(~arid, nrow = 1) + ylab('Percent of Peaks') + xlab('') 
all_near_plot

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

# make he grange objects
arid1a_gr <- with(arid1a_peaks, GRanges(seqnames = seqnames, IRanges(start, end), strand = '*', name = name))
arid1b_gr <- with(arid1b_peaks, GRanges(seqnames = seqnames, IRanges(start, end), strand = '*', name = name)) 
arid2_gr  <- with(arid2_peaks, GRanges(seqnames = seqnames, IRanges(start, end), strand = '*', name = name) ) 


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

pick_random_hoods <- function(gr, expr, p_expr) {
  # need a way to get random expression for a set of peaks (arid1a, arid1b, arid2) .
  # this function gets called during compare_expression_in_hood and requires the 
  # output of neighborhood_expreesion to give the distirbution of number of gnees
  # I'll randomly select that number of genes from the expr data and return the data
  gene_num_dist <- unlist(lapply(p_expr$genes_in_hood, function(x) length(unlist(x) ) )) 
  num_to_sample <- sample(gene_num_dist, length(gr), replace =T )
  # get a list of genes corresponding to these numbers
  vals <- lapply(num_to_sample, function(x) mean(expr[sample(1:nrow(expr), size = x, replace =F), 'log2FoldChange'] ))
  gr %<>% as.data.frame() %>% tbl_df() %>%
    mutate(mean_hood_expr = unlist(vals), genes_in_hood = rep('NA')) %>%
    select(name, genes_in_hood, mean_hood_expr) 
  return(gr) 
}

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
 
all_distances %>% 
 mutate(distance = distance/1000) %>%
 ggplot(aes(x = as.factor(distance), y = abs(mean_hood_expr), fill = random, color = random) ) + 
  geom_boxplot(notch = T, width = 0.9, color = 'grey20', position = position_dodge(0.9), outlier.shape = NA, outlier.size = NA) + 
  facet_wrap(~arid, nrow =1, scale = 'free') + theme_paper() + scale_fill_manual(values = c('steelblue', 'grey30') ) + 
  scale_y_continuous(limits = c(0, 1.5)) + 
  theme(axis.text.x = element_text(color = 'grey10'), 
        axis.ticks = element_blank(), 
        panel.background = element_rect(fill = NA, color = 'grey30')) +
  xlab('Distance Surrounding Peak (kb)') + ylab('Abs Mean Expression\nof Genes in Neighborhood')

all_distances %>%
  group_by(arid, distance) %>%
  do(tidy(wilcox.test(abs(mean_hood_expr) ~ random, data = .) ))


arid1a_plots <-  arid1a_peaks_plus %>% ggplot(aes(x = svm, y = mean_hood_expr )) + geom_jitter(alpha = 0.1) + geom_boxplot(fill =NA)
arid1b_plots <-  arid1b_peaks_plus %>% ggplot(aes(x = svm, y = mean_hood_expr)) + geom_jitter(alpha = 0.1) + geom_boxplot(fill = NA) 
arid2_plots <- arid2_peaks_plus %>% ggplot(aes(x = svm, y = mean_hood_expr)) + geom_jitter(alpha = 0.1) + geom_boxplot(fill =NA)
grid.arrange(arid1a_plots, arid1b_plots, arid2_plots, ncol =3) 

#compare numbers wacross groups. 
# may want to function this up so that I can do the different values of 'neighborhood'
#stogram(aes(y=..count../sum(..count..)), binwidth=0.1, alpha  = 0.25)  + 
all_arids_forhist <- rbind(arid1a_peaks_plus, arid1b_peaks_plus, arid2_peaks_plus) 
aridhist <- all_arids_forhist %>% ggplot(aes(x=mean_hood_expr, color = svm, fill = svm) ) +
  geom_density(alpha = 0.25) + 
  scale_fill_manual(values = c('grey30', 'steelblue') ) + 
  scale_color_manual(values = c('grey30', 'steelblue')) + 
  theme_paper() + theme(panel.grid.major.x = element_blank(), 
                        axis.ticks.y = element_blank(), 
                        panel.grid.major.y = element_line(color = 'grey70') ) + 
  xlab('') + ylab('') +facet_wrap(~arid, ncol =1) 
aridhist

wilcox.test(abs(arid1a_peaks[arid1a_peaks$svm == 'Multiple', 'mean_hood_expr']), abs(arid1a_peaks[arid1a_peaks$svm == 'Single', 'mean_hood_expr']) )
wilcox.test(abs(arid1b_peaks[arid1b_peaks$svm == 'Multiple', 'mean_hood_expr']), abs(arid1b_peaks[arid1b_peaks$svm == 'Single', 'mean_hood_expr']) )
wilcox.test(abs(arid2_peaks[arid2_peaks$svm == 'Multiple', 'mean_hood_expr']), abs(arid2_peaks[arid2_peaks$svm == 'Single', 'mean_hood_expr']) ) 

g1 <- ggplot(arid2_peaks, aes(x= log10(distance), y = abs(mean_hood_expr), color = svm)) + geom_point(alpha = 0.3)  +
  scale_color_manual(values = c('red2', 'steelblue')) + 
  theme(plot.margin=unit(c(0,0,1,1), 'lines') )
g2 <- ggplot(arid2_peaks, aes(x = log10(distance), color = svm)) + geom_density() + 
  scale_color_manual(values = c('red2', 'steelblue')) + xlab('') + 
  theme(legend.position = 'none', axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(1,1,-1,1), 'lines') ) 
g3 <- ggplot(arid2_peaks, aes(x = abs(mean_hood_expr), color = svm)) + geom_density() + coord_flip() + 
  theme(legend.position = 'none', axis.title.y = element_blank(), axis.text.y = element_blank(), 
        axis.ticks = element_blank(), plot.margin = unit(c(1, 1, 1, -1), 'lines') ) + 
  scale_color_manual(values = c('red2', 'steelblue'))

g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

legend2<- g_legend(g1 + theme(legend.key.size = unit(1, 'cm') ) )

grid.arrange(g2, legend2, g1+theme(legend.position = 'none'), g3, widths = c(0.8, 0.2), heights = c(0.2, 0.8) )
empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
       panel.background=element_blank(), 
       axis.text=element_blank(),           
       axis.title.x=element_blank(), 
       panel.grid = element_blank() )

# might be better to shuffle the genes around - 
# selecting a random number of genes from gencode for each group and giving expression values for those

arid1b_peaks_noarid2change <- arid1b_peaks[!arid1b_peaks$nearestTSS %in% arid2_changed$rowname, ]
arid1b_peaks_noarid2change <- arid1b_peaks_noarid2change[arid1b_peaks_noarid2change$nearestTSS %in% arid1b_changed$rowname, ] 

arid1b_peaks_noarid2change %>% 
  group_by(type, svm) %>% 
  summarise(count = n() )

# to the arid peaks list add if nearestTSS is significantly changed
# No obvious diff between single vs multiple bound peaks as far as 
# genomic distribution of peaks near significantly changed genes. 
isSigFun <- function(peaks, expr) { 
  peaks %>% 
    mutate(isSig = ifelse(nearestTSS %in% expr$rowname, 1, 0) ) %>%
    group_by(type, svm, isSig) %>% 
    summarise(numsig = n() ) %>%
    mutate(percent = numsig/sum(numsig) )
  }
# Is there evidence that single vs multiple has more signficant changed expression 
# - hypothesis would be no, we see  more evidence for arid1b and arid2 being necessary for repression - and in competition arid1a or arid2 loss has an effect 
# - therefore you might expect if anything multiple to be enriched. Formally show this. 

# compare percent of signficant genes changed in each of the two categoreis. 

# Figure 8-1 C
arid1a_is_sig <- arid1a_peaks %>% isSigFun(arid1a_changed) %>% mutate(arid = 'Arid1a') 
arid1b_is_sig <- arid1b_peaks %>% isSigFun(arid1b_changed) %>% mutate(arid = 'Arid1b') 
arid2_is_sig <- arid2_peaks %>% isSigFun(arid2_changed ) %>% mutate(arid = 'Arid2') 
arid_sig_comb <- rbind(arid1a_is_sig, arid1b_is_sig, arid2_is_sig) 
arid_sig_comb %>% 
  filter(isSig == 1) %>% 
  ggplot(aes(x = svm, y = numsig, fill = arid)) +
  geom_bar(position = 'dodge', stat = 'identity') + 
  geom_bar(position = 'dodge', stat = 'identity', color = 'grey30', show_guide =F) +
  facet_wrap(~type, nrow =1) +
  theme_paper()  +
  theme(panel.background = element_rect(color = 'grey30', fill = NA), 
        panel.grid = element_blank(), 
        legend.title = element_blank() ) + scale_y_continuous(expand = c(0,0) ) + 
  ylab('Number of Peaks Nearest \nto Significantly Changed Genes') + xlab('') + 
  scale_fill_manual(values = cols) +xlab('') 

arid1b_peaks %>% mutate(isSig = ifelse(nearestTSS %in% arid1b_changed$rowname, 1, 0) ) %>%
  left_join(arid1b_all, by = c('nearestTSS' = 'rowname') )  %>% 
  filter(padj < 0.05) %>% 
  mutate(direction = ifelse(log2FoldChange < 0, 'Down', 'Up') ) %>%
  group_by(type, svm, direction) %>% 
  summarise(count = n() ) %>% 
  mutate(percent  = count/sum(count) ) %>% 
  ggplot(aes(x = type, y = percent, fill = direction)) + geom_bar(position = 'fill', stat = 'identity') + theme_paper() + 
    facet_wrap(~svm, ncol = 1) 

#when Arid1b and Arid2 co-repress a gene directly, where are they bound (together or separately) 
#Filter the peakl ists by the genes repressed by both
rep1b <- arid1b_changed %>% filter(log2FoldChange > 0) %>% select(rowname, log2FoldChange, padj, baseMean) 
rep2 <- arid2_changed %>% filter(log2FoldChange > 0) %>% select(rowname, log2FoldChange, padj, baseMean) 
repressed_in_both <- inner_join(rep1b, rep2, by = 'rowname', x_suffix ='Arid1b', y_suffix = 'Arid2')  
act1b <- arid1b_changed %>% filter(log2FoldChange < 0) %>% select(rowname, log2FoldChange, padj, baseMean) 
act2  <- arid2_changed %>% filter(log2FoldChange < 0 ) %>% select(rowname, log2FoldChange, padj, baseMean)
act_in_both <- inner_join(act1b, act2, by = 'rowname', x_suffix = 'Arid1b', y_suffix = 'Arid2' )

arid1b_peaks_repr_direct <- arid1b_peaks %>% 
  filter(nearestTSS %in% repressed_in_both$rowname) %>% 
  mutate(arid = 'Arid1b') 
  
arid2_peaks_repr_direct <- arid2_peaks %>% 
  filter(nearestTSS %in% repressed_in_both$rowname) %>%
  mutate(arid ='Arid2') 

arid1b_peaks_act_direct <- arid1b_peaks %>%
  filter(nearestTSS %in% act_in_both$rowname) %>%
  mutate(arid = 'Arid1b') 

arid2_peaks_act_direct <- arid2_peaks %>%
  filter(nearestTSS %in% act_in_both$rowname) %>%
  mutate(arid = 'Arid2') 

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

# Do those peaks overlap and where are they. 

#full tables 
arid1b_drb_ft <- arid1b_peaks_repr_direct %>% 
  filter(nearestTSS %in% direct_repr_both) 
# 233 peaks covering those 114 genes

arid2_drb_ft <- arid2_peaks_repr_direct %>% 
  filter(nearestTSS %in% direct_repr_both)
# 314 peaks coverign these 114 genes

arid1b_dab_ft <- arid1b_peaks_act_direct %>%
  filter(nearestTSS %in% direct_act_both)
# 72 peaks cover the 44 genes

arid2_dab_ft <- arid2_peaks_act_direct %>%
  filter(nearestTSS %in% direct_act_both)
# 82 peaks cover the 44 genes


arid1b_drb_summ <- arid1b_drb_ft %>% 
  group_by(type, svm ) %>% 
  summarise(count = n() ) %>%
  mutate(arid = 'Arid1b', dir = 'Repressed') 

arid2_drb_summ <- arid2_drb_ft %>% 
  group_by(type, svm) %>%
  summarise(count = n() ) %>% 
  mutate(arid = 'Arid2', dir = 'Repressed') 

arid1b_dab_summ <- arid1b_dab_ft %>% 
  group_by(type, svm ) %>% 
  summarise(count = n() ) %>%
  mutate(arid = 'Arid1b', dir = 'Activated') 

arid2_dab_summ <- arid2_dab_ft %>% 
  group_by(type, svm) %>%
  summarise(count = n() ) %>% 
  mutate(arid = 'Arid2', dir = 'Activated') 

# This doens't really answer anything. 
# distal sites are 50/50
# promoters are mayb 60/40 (multipole) 
# other two are more single
# no real diff between 1b and 2

combined_arid1b_2_summ <- rbind(arid1b_drb_summ, arid2_drb_summ, arid1b_dab_summ, arid2_dab_summ)  %>%
  ungroup() %>%
  arrange(desc(type), desc(count)) %>%
  mutate(type = factor(type, levels = c('Promoter', 'Exon', 'Distal', 'Intron') ) )
combined_arid1b_2_summ %>% 
  ggplot(aes(x=svm, y = count, fill = dir)) + geom_bar(position = 'dodge', stat = 'identity') +
  facet_grid(arid~type)   +theme_paper() + 
  scale_fill_manual(values = c('red3', 'steelblue') ) 

#same as above but put the fill as the type rather than direction 
combined_arid1b_2_summ %>%
  ggplot(aes(x=dir, y = count, fill = type, order = type)) + geom_bar(position = 'fill', stat = 'identity') + 
  geom_bar(position = 'fill', stat = 'identity', color ='grey20', show_guide = F) +
  facet_grid(arid~svm) + theme_paper() + scale_fill_manual(values = pal1) +
  theme(axis.ticks = element_blank(), 
        axis.text = element_text(color = 'grey20') ) +
  xlab('') + ylab('Percent of Peaks') 



#################################################################################################
#what are the peaks that overlap. 
# basically want to know if when regulating these genes they are bound at a particular type of peak
# need to first collect the peaks that overlap and are the direct regulating peaks. 
# create GRange objects for each of these
arid1b_drb_peaks <- with(arid1b_drb_ft, GRanges(seqnames = seqnames, IRanges(start = start, end = end), strand = '*') )
elementMetadata(arid1b_drb_peaks) <- arid1b_drb_ft[,4:ncol(arid1b_drb_ft)]
arid2_drb_peaks <- with(arid2_drb_ft, GRanges(seqnames = seqnames, IRanges(start, end), strand = '*') ) 
elementMetadata(arid2_drb_peaks) <- arid2_drb_ft[,4:ncol(arid2_drb_ft) ] 
arid1b_dab_peaks <- with(arid1b_dab_ft, GRanges(seqnames, IRanges(start, end), strand = '*'))
elementMetadata(arid1b_dab_peaks) <- arid1b_dab_ft[,4:ncol(arid1b_dab_ft)] 
arid2_dab_peaks <- with(arid2_dab_ft, GRanges(seqnames, IRanges(start, end), strand = '*') )
elementMetadata(arid2_dab_peaks) <- arid2_dab_ft[,4:ncol(arid2_dab_ft) ] 

#compare the overlap of the two peak sets
ovp_rep <- intersect(arid1b_drb_peaks, arid2_drb_peaks) 
ovp_rep$name <- unlist(lapply(1:length(ovp_rep), function(x) paste('ovpeak_', x, sep ='') ))
ovp_act <- intersect(arid1b_dab_peaks, arid2_dab_peaks) 
ovp_act$name <- unlist(lapply(1:length(ovp_act), function(x) paste('ovpeak_', x, sep ='') ) )

# 82 of the peaks overlap. 
ovp_rep_anno <- annotatePeaks(peaks = ovp_rep, gencode, window =2000)
ovp_rep_anno$group <- 'Repressed'
ovp_act_anno <- annotatePeaks(peaks = ovp_act, gencode, window = 2000) 
ovp_act_anno$group <- 'Activated' 
overlapped_annotated <- c(ovp_rep_anno, ovp_act_anno) 

tmp <- overlapped_annotated %>% as.data.frame() %>% tbl_df() %>% 
  group_by(type, group) %>% 
  summarise(count = n() ) 

tmp %>% 
  spread(group, count) %>%
  select(-type) %>%
  fisher.test(.) 

# no difference between groups when considering the 82 peaks that are overlapping between Arid1b and Arid2 and regulate genes that are nearest to these peaks in both sets. 
# this is the most wittled down set of genes. 

tmp %>%
  ggplot(aes(x = group, y = count, fill = type)) + geom_bar(stat = 'identity', position = 'fill') +
    geom_bar(stat = 'identity', position = 'fill', color = 'grey20', show_guide = F) + 
    theme_paper() + theme(axis.ticks = element_blank(), axis.text = element_text(color = 'grey20') )+ 
    scale_fill_brewer(palette = 'Set2') + 
    xlab('') + ylab('Fraction of Peaks') 
  




# The derepressed targets with overlapping peaks are mostly promoter
# shoudl compare to upregulated directly co=ocupied

# how bout peaks that are not overlapping
# Takea  step back and look at all peaks that are associated with double repressed or doubly activated genes. 
# This was used for an older version of 8C
# difference between it and below is this is a more relaxed(and larger) peak set
# that required the gene be regulated by both, but only one must be assigned to it
# 
arid1b_pas <- arid1b_peaks_act_direct %>% 
  group_by(type) %>% 
  summarise(count = n() ) %>% 
  mutate(group = 'Arid1b', dir = 'Activated') 
arid2_pas <- arid2_peaks_act_direct %>% 
  group_by(type) %>% 
  summarise(count = n() ) %>%
  mutate(group = 'Arid2', dir = 'Activated') 
arid1b_prs <- arid1b_peaks_repr_direct %>% 
  group_by(type) %>% 
  summarise(count = n() ) %>% 
  mutate(group = 'Arid1b', dir = 'Repressed') 
arid2_prs <- arid2_peaks_repr_direct %>% 
  group_by(type) %>% 
  summarise(count = n() ) %>%
  mutate(group = 'Arid2', dir = 'Repressed') 
all_direct_summ <- rbind(arid1b_pas, arid2_pas, arid1b_prs, arid2_prs) %>%
  mutate(type = factor(type, levels = c('Promoter', 'Exon', 'Distal', 'Intron') ) )
all_direct_summ %>% 
  ggplot(aes(x = dir, y = count, fill = type, order = type)) + 
    geom_bar(position = 'fill', stat = 'identity') +
    geom_bar(position = 'fill', stat = 'identity', color = 'grey20', show_guide=F)+ 
    facet_wrap(~group,  nrow = 1)  + 
    theme_paper() + 
    theme(axis.text = element_text(color = 'grey20'), 
          axis.ticks = element_blank() ) + 
    scale_fill_manual(values = pal1) + scale_y_continuous(expand = c(0,0) ) + 
    xlab('') + ylab('Proportion of peaks')
    
all_direct_summ %>% 
  filter(group == 'Arid1b') %>%
  spread(dir, count)  %>%  
  select(-type, -group) %>%  
  do(tidy(chisq.test(cbind(.$Activated, .$Repressed)) ) ) 

all_direct_summ %>% 
  filter(group == 'Arid2') %>% 
  spread(dir, count) %>% print() %>% 
  select(-type, -group) %>% 
  do(tidy(chisq.test(cbind(.$Activated, .$Repressed)))) 


#j peaks directly acted on by both. 
# this is what I actually used for figure 8C
both_act_direct <- intersect( both_act, intersect(arid1b_peaks$nearestTSS, arid2_peaks$nearestTSS) )
# 44 genes associated with this
arid1b_dbl_dir_act <- arid1b_act_direct %>% 
  filter(nearestTSS %in% both_act_direct) 
arid2_dbl_dir_act <- arid2_act_direct %>% 
  filter(nearestTSS %in% both_act_direct) 
# 72 peaks
both_repr_direct <- intersect(both_repr, intersect(arid1b_peaks$nearestTSS, arid2_peaks$nearestTSS))
arid1b_dbl_dir_rep <- arid1b_repr_direct %>% 
  filter(nearestTSS %in% both_repr_direct) 
arid2_dbl_dir_rep <- arid2_repr_direct %>% 
  filter(nearestTSS %in% both_repr_direct) 


arid1b_pas <- arid1b_dbl_dir_act%>% 
  group_by(type) %>% 
  summarise(count = n() ) %>% 
  mutate(group = 'Arid1b', dir = 'Activated') 
arid2_pas <- arid2_dbl_dir_act%>% 
  group_by(type) %>% 
  summarise(count = n() ) %>%
  mutate(group = 'Arid2', dir = 'Activated') 
arid1b_prs <- arid1b_dbl_dir_rep %>% 
  group_by(type) %>% 
  summarise(count = n() ) %>% 
  mutate(group = 'Arid1b', dir = 'Repressed') 
arid2_prs <- arid2_dbl_dir_rep %>% 
  group_by(type) %>% 
  summarise(count = n() ) %>%
  mutate(group = 'Arid2', dir = 'Repressed') 

all_direct_summ <- rbind(arid1b_pas, arid2_pas, arid1b_prs, arid2_prs) %>%
  mutate(type = factor(type, levels = c('Promoter', 'Exon', 'Distal', 'Intron') ) )
all_direct_summ %>% 
  ggplot(aes(x = dir, y = count, fill = type, order = type)) + 
  geom_bar(position = 'fill', stat = 'identity') +
  geom_bar(position = 'fill', stat = 'identity', color = 'grey20', show_guide=F)+ 
  facet_wrap(~group,  nrow = 1)  + 
  theme_paper() + 
  theme(axis.text = element_text(color = 'grey20'), 
        axis.ticks = element_blank() ) + 
  scale_fill_manual(values = pal1) + scale_y_continuous(expand = c(0,0) ) + 
  xlab('') + ylab('Proportion of peaks')

all_direct_summ %>% 
  filter(group == 'Arid1b') %>%
  spread(dir, count)  %>%  
  select(-type, -group) %>%  
  do(tidy(chisq.test(cbind(.$Activated, .$Repressed)) ) ) 

all_direct_summ %>% 
  filter(group == 'Arid2') %>% 
  spread(dir, count) %>% print() %>% 
  select(-type, -group) %>% 
  do(tidy(chisq.test(cbind(.$Activated, .$Repressed)))) 

