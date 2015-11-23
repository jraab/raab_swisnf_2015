# Script to create Figure 4 C, D, E
# -----------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(magrittr)
library(GenomicRanges)
library(VennDiagram)
library(Hmisc)
library(rhdf5)
library(gridExtra)
library(broom)
source('code/function_defs/chipseq_fig_fun.R') # for bed2gr
source('code/util/theme_paper.R') 

# Important Directories
# -----------------------------------------------------------------------------
output_figs <- 'figures/'
chip_peaks <- 'output/macs_peaks/cleaned/'
# Requires files from output of deeptools ARID analysis to be present in 
# output/plots/dt_matfiles/
# output/encode_coverages/coverages - must have h3k27ac, h3k4me3, h3k4me1 hdf5 files

# Load Data 
# ----------------------------------------------------------------------------
# files for comparing Arids
arid1a <- read.csv('output/macs_peaks/cleaned/arid1a_annotated.csv')
arid1b <- read.csv('output/macs_peaks/cleaned/arid1b_annotated.csv')
arid2  <- read.csv('output/macs_peaks/cleaned/arid2_annotated.csv')
snf5   <- read.csv('output/macs_peaks/cleaned/snf5_annotated.csv') 
all <- rbind(arid1a, arid1b, arid2, snf5) 
all %<>% separate(name, into =c('Arid', 'Name'), sep = '_', extra = 'merge') 
all$type <- factor(as.character(all$type), levels = c('Promoter', 'Exon', 'Intron', 'Distal') )
all$Arid <- capitalize(all$Arid) 
all %<>% arrange(type) 
# objects for doing Venn Overlaps
snf5_gr <- 'output/macs_peaks/cleaned/snf5_peaks_cf.bed'
arid1a_gr <- bed2gr(paste0(chip_peaks, 'arid1a_peaks_cf.bed'), stranded = F )
arid1b_gr <- bed2gr(paste0(chip_peaks, 'arid1b_peaks_cf.bed'), stranded = F )
arid2_gr  <- bed2gr(paste0(chip_peaks, 'arid2_peaks_cf.bed' ), stranded = F )
combined_arids <- union(arid1a_gr, union(arid1b_gr, arid2_gr) ) 
snf5_gr   <- bed2gr(paste0(chip_peaks, 'snf5_peaks_cf.bed'  ), stranded = F)
snyder_snf5 <- 'data/external/snyderData/wgEncodeSydhTfbsHelas3Ini1IggmusPk.narrowPeak.gz'
snyder_gr <- bed2gr(snyder_snf5, stranded = F)  



# Genomic Location- Figure 4C
# ----------------------------------------------------------------------------
pal1 <- rev(c('#ef3b2c', '#fc9272', '#2171b5', '#6baed6'))
pal2 <- rev(c('#238b45', '#74c476', '#2171b5', '#6baed6'))

summ <- all %>% group_by(Arid, type) %>% summarise(number= n() ) %>% mutate(percent = number/sum(number) *100  )
g <- ggplot(summ, aes(x = toupper(Arid), y = percent, fill = type) ) + geom_bar(stat = 'identity') 
g <- g + geom_bar(stat = 'identity', color = 'grey30', show_guide = F) 
g <- g+ theme_figure() 
g <- g + scale_fill_manual(values = pal1) +xlab('') + ylab('Percent') 
g <- g + scale_y_continuous(expand = c(0,0) )
g <- g + theme(axis.text = element_text(color = 'grey20', size = 18),
               axis.ticks = element_blank())
g

#save figure 
ggsave(paste0(output_figs, 'fig3c_genomic_distribution.pdf'), g) 

# Venn Diagrams - Supplemental Figure 4 -1
# ---------------------------------------------------------------------------
# 2 way with snf 5
# red blue green purple grey
venn <- c('#e34a33','#2b8cbe', '#2ca25f', '#756bb1', '#636363')
pdf('figures/supplemental/3A_arid1a_snf.pdf')
draw2way(arid1a_gr, snf5_gr, c('ARID1A', 'SNF5'), colors = c(venn[1], venn[4]) )
dev.off()
pdf('figures/supplemental/3B_arid1b_snf.pdf')
draw2way(arid1b_gr, snf5_gr, c('ARID1B', 'SNF5'), colors = c(venn[2], venn[4]) )
dev.off()
pdf('figures/supplemental/3C_arid2_snf.pdf') 
draw2way(arid2_gr, snf5_gr, c('ARID2', 'SNF5'), colors = c(venn[3], venn[4]) )
dev.off()
pdf('figures/supplemental/3D_arids_snf.pdf')
draw2way(combined_arids, snf5_gr, c('All ARIDs', 'SNF5'), colors= c(venn[5],venn[4]))
dev.off()
pdf('figures/supplemental/3E_snf5s.pdf')
draw2way(snyder_gr, snf5_gr, c('Euskirchen', 'This Study'), colors = c(venn[5], venn[4]) ) 
dev.off()

# Figures 4D,E 
# --------------------------------------------------------------------------
# This relies on the signal calculated by deeptools 
# It is a quantification of values heatmape in 4A,D


# Function to read in output from deeptools matrix and split 
# into two groups based on the bed file use for deep tools 
# return a data frame with the rowmean of the signal
split_vals <- function(filename) { 
  val_df <- read.table(filename, skip = 1)
  val_df <- val_df[,8:ncol(val_df)]
  f <- readLines(filename) 
  split_rows <- which(grepl(pattern = '#', f) )
  active <- val_df[1:split_rows[1]-1, ]
  poised <- val_df[(split_rows[1]+1):nrow(val_df), ]
  # check that I Only have two groups
  if (length(split_rows) != 2 ) { 
    print('ERROR: Should only be two groups to split') 
  }
  df <- data.frame(Vals = c(rowMeans(active), rowMeans(poised)), 
                   Group = c(rep('Active', nrow(active)), rep('Poised', nrow(poised) ) ) )
  return(df) 
  
}

# Slight change from above to look at active vs inactive promoters
split_vals_genes <- function(filename) { 
  val_df <- read.table(filename, skip = 1)
  val_df <- val_df[,8:ncol(val_df)]
  f <- readLines(filename) 
  split_rows <- which(grepl(pattern = '#', f) )
  active <- val_df[1:split_rows[1]-1, ]
  poised <- val_df[(split_rows[1]+1):nrow(val_df), ]
  # check that I Only have two groups
  if (length(split_rows) != 2 ) { 
    print('ERROR: Should only be two groups to split') 
  }
  df <- data.frame(Vals = c(rowMeans(active), rowMeans(poised)), 
                   Group = c(rep('Active', nrow(active)), rep('Inactive', nrow(poised) ) ) )
  return(df) 
}

# Do the work
# arid enhancer files (4E)
arid1a <- 'output/dt_matfiles/enhancers_arid1a_mat'
arid1b <- 'output/dt_matfiles/enhancers_arid1b_mat'
arid2 <- 'output/dt_matfiles/enhancers_arid2_mat'
snf5 <- 'output/dt_matfiles/enhancers_snf5_mat'



arid1a_df <- split_vals(arid1a) 
arid1b_df <- split_vals(arid1b)
arid2_df <- split_vals(arid2) 
snf5_df  <- split_vals(snf5) 
arid1a_df$name <- rep('ARID1A')
arid1b_df$name <- rep('ARID1B') 
arid2_df$name <- rep('ARID2') 
snf5_df$name  <- rep('SNF5') 
all <- rbind(arid1a_df, arid1b_df, arid2_df, snf5_df) 
g <- ggplot(all, aes(x=Group, y = log2(Vals)) ) + geom_boxplot(notch =T, fill = 'grey60', width = 0.5) 
g <- g + facet_wrap(~name, ncol = 4) + theme_figure()
g <- g + xlab('') + ylab('Log2 Arid Signal') 
g <- g + theme(axis.text = element_text(color = 'grey20', angle = 30, hjust = 1, vjust =1) ) 
g
ggsave('figures/fig3e_aridsig_act_v_pois_enh.pdf', g)

# Is differnce signfiicant - yes
arid1a_df %>% group_by(name) %>%
  wilcox.test(Vals~Group, data = .) %>% 
  tidy()

arid1b_df %>% group_by(name) %>% 
  wilcox.test(Vals~Group, data = .) %>%
  tidy()

arid2_df %>% group_by(name) %>% 
  wilcox.test(Vals~Group, data = .) %>% 
  tidy()

# Do test for promoters ( 4D ) 
arid1a_g <- 'output/dt_matfiles/genes_arid1a_mat'
arid1b_g <- 'output/dt_matfiles/genes_arid1b_mat'
arid2_g  <- 'output/dt_matfiles/genes_arid2_mat'
snf5_g   <- 'output/dt_matfiles/genes_snf5_mat'

arid1a_gdf <- split_vals_genes(arid1a_g) 
arid1b_gdf <- split_vals_genes(arid1b_g)
arid2_gdf <- split_vals_genes(arid2_g) 
snf5_gdf  <- split_vals_genes(snf5_g) 
arid1a_gdf$name <- rep('ARID1A')
arid1b_gdf$name <- rep('ARID1B') 
arid2_gdf$name <- rep('ARID2') 
snf5_gdf$name   <- rep('SNF5') 
allg <- rbind(arid1a_gdf, arid1b_gdf, arid2_gdf, snf5_gdf)  
h <- ggplot(allg, aes(x=Group, y = log2(Vals+1)) ) + geom_boxplot(notch =T, width = 0.5, fill = 'grey60') 
h <- h + facet_wrap(~name, ncol = 4) + theme_figure()
h <- h + xlab('') + ylab('Log2 Arid Signal') 
h <- h + theme(axis.text = element_text(color = 'grey20', angle = 30, hjust = 1, vjust =1) ) 
h

ggsave('figures/fig3d_aridsig_act_v_pois_gene.pdf', h)

arid1a_gdf %>% group_by(name) %>%
  wilcox.test(Vals~Group, data = .) %>% 
  tidy()

arid1b_gdf %>% group_by(name) %>% 
  wilcox.test(Vals~Group, data = .) %>%
  tidy()

arid2_gdf %>% group_by(name) %>% 
  wilcox.test(Vals~Group, data = .) %>% 
  tidy()

# Supplement 4-2
arid1a <- read.csv('output/macs_peaks/cleaned/arid1a_annotated.csv') %>% 
  mutate(arid = 'ARID1A' ) %>% filter(svm == 'Single') %>% select(state, arid)
arid1b <- read.csv('output/macs_peaks/cleaned/arid1b_annotated.csv') %>% 
  mutate(arid = 'ARID1B' ) %>% filter(svm == 'Single') %>% select(state, arid)
arid2  <- read.csv('output/macs_peaks/cleaned/arid2_annotated.csv') %>%
  mutate( arid = 'ARID2' ) %>% filter(svm == 'Single') %>% select(state, arid)
multiple <- read.csv('output/macs_peaks/cleaned/multi_bound_peaks.csv', header =F) %>% 
  mutate(arid = 'Multiple') 
colnames(multiple) <- c('chr', 'start', 'end', 'name', 'nearestTSS', 'dist', 'state', 'arid') 
multiple <- multiple %>% select(state, arid)

all <- rbind(arid1a, arid1b, arid2, multiple) 
all <- all %>% 
  group_by(arid, state) %>% 
  summarise(count = n() ) %>% 
  mutate( percent = count/sum(count)  )  %>% 
  filter(!state %in% 'NONE')

a <- all %>% 
  separate(state, into = c('number', 'name'), sep ='_', extra = 'merge' ) %>% 
  arrange(as.numeric(number)) %>% mutate(newnames=paste(number, name, sep ='_') )

lvls <- rev(as.character(unique(a$newnames) ) )

p <- a %>%
  mutate(newnames2 = factor(newnames, levels = lvls) )%>% 
  ggplot(aes(x = arid, y = newnames2, fill = percent)) + 
  geom_tile(color = 'grey30') + scale_fill_gradient(low = 'white', high = 'red2')  + 
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 18) ) + xlab('') + ylab('') 
p
ggsave('figures/supplemental/s3_2_chromstate.pdf')

# Supplement 4-3
# ----------------------------------------------------------------------------
stubs <- c('arid1a', 'arid1b', 'arid2') 
ends <- c('_p_H3k04me1StdAln_s_coverage.h5', 
          '_p_H3k4me3StdAln_s_coverage.h5', 
          '_p_H3k27acStdAln_s_coverage.h5')
# H3K4me1
me1 <- plot_arids(ends[1], stubs) + theme(legend.position = 'none') 
# H3K4me3
me3 <- plot_arids(ends[2], stubs)
# H3K27ac
k27 <- plot_arids(ends[3], stubs) + theme(legend.position = 'none') 
pdf('figures/supplemental/fig4_3_metaplots.pdf', height = 4, width = 12)
grid.arrange(me1, me3+xlab('Distance from Peak Center (Kb)'), k27, nrow =1)
dev.off()
