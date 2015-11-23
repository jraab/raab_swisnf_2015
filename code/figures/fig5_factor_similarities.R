# Figure 6 
# -----------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(magrittr)
library(GenomicRanges)
library(Hmisc)
library(data.table)
library(RColorBrewer)
library(gplots)
library(rhdf5) 
source('code/function_defs/fig5_helpers_func.R')
source('code/function_defs/cleanNames.R')
source('code/util/theme_paper.R')

# Load Data - Information about alonve vs multi peaks
# ------------------------------------------------------------------------------
# Peak information
arid1a <- read.table('output/macs_peaks/cleaned/arid1a_annotated.csv', sep =',', header =T) 
arid1a_single <- arid1a[arid1a$svm == 'Single', ]
arid1b <- read.table('output/macs_peaks/cleaned/arid1b_annotated.csv', sep =',', header =T)  
arid1b_single <- arid1b[arid1b$svm == 'Single', ]
arid2 <- read.table('output/macs_peaks/cleaned/arid2_annotated.csv', sep = ',', header =T) 
arid2_single <- arid2[arid2$svm == 'Single', ]
multi <- read.table('output/macs_peaks/cleaned/multi_bound_peaks.csv', sep = ',', header =F) # multi doesn't seem to have a name
colnames(multi) <- c('seqnames', 'start', 'end','name')  

# Enrichment Information # These are for all peaks _ so neeed filtered by 'Single' peaks
arid1a_enrich <- read.csv('output/encode_coverages/Arid1a_enrichments_bypeaks.csv')
arid1b_enrich <- read.csv('output/encode_coverages/Arid1b_enrichments_bypeaks.csv') 
arid2_enrich  <- read.csv('output/encode_coverages/Arid2_enrichments_bypeaks.csv') 
multi_enrich  <- read.csv('output/encode_coverages/Multi_enrichments_bypeaks.csv') %>%
  mutate(peakname = 'Multi')

arid1a_enrich <- arid1a_enrich[arid1a_enrich$gene %in% arid1a_single$name, ] %>% 
  mutate(peakname = 'ARID1A')
arid1b_enrich <- arid1b_enrich[arid1b_enrich$gene %in% arid1b_single$name, ] %>% 
  mutate(peakname = 'ARID1B')
arid2_enrich  <- arid2_enrich[arid2_enrich$gene %in% arid2_single$name, ] %>% 
  mutate(peakname = 'ARID2')


# Correlation Information 
arid1a_corr <- read.table('output/encode_coverages/arid1a_factors_cor.tab', header =T, sep = '\t')
arid1b_corr <- read.table('output/encode_coverages/arid1b_factors_cor.tab', header =T, sep = '\t') 
arid2_corr  <- read.table('output/encode_coverages/arid2_factors_cor.tab', header = T, sep ='\t') 
multi_corr <- read.table('output/encode_coverages/multi_factors_cor.tab', header =T, sep = '\t') 

# ------------------------------------------------------------------------------
# Combine and summarise the enrichment values 
full = rbind(arid1a_enrich, arid1b_enrich, arid2_enrich, multi_enrich)
full_summary <- full %>% 
  gather(signame , enrichment, -gene, -peakname) %>% 
  group_by(peakname, signame) %>%
  summarise(avg=mean(enrichment) )

arid1a_sum <- full_summary[full_summary$peakname == 'ARID1A', ]
arid1b_sum <- full_summary[full_summary$peakname == 'ARID1B', ]
arid2_sum <- full_summary[full_summary$peakname == 'ARID2',] 
multi_sum <- full_summary[full_summary$peakname == 'Multi', ]

# This saves the plots for Supplemental Figure 6-1
# ------------------------------------------------------------------------------
save_info(arid1a_corr, arid1a_sum, 'arid1a') 
save_info(arid1b_corr, arid1b_sum, 'arid1b') 
save_info(arid2_corr, arid2_sum, 'arid2') 
save_info(multi_corr, multi_sum, 'multi')

## Uncomment this section if you wish to run interactively and view plots
#save_info(arid1a_corr, arid1a_sum, 'arid1a', print =TRUE) 
#save_info(arid1b_corr, arid1b_sum, 'arid1b', print = TRUE) 
#save_info(arid2_corr, arid2_sum, 'arid2', print = T) 
#save_info(multi_corr, multi_sum, 'multi', print =T) 

# Network Plot Showing shared factors between clusters of different peak types
# -----------------------------------------------------------------------------
# Create a network plot showing relationship between clusters. 
#Share category
arid1a_clus <- read.csv('output/supplemental_data/fig6_arid1a_clusterassign.csv') 
arid1b_clus <- read.csv('output/supplemental_data/fig6_arid1b_clusterassign.csv') 
arid2_clus <- read.csv('output/supplemental_data/fig6_arid2_clusterassign.csv') 
multi <- read.csv('output/supplemental_data/fig6_multi_clusterassign.csv') 

all <- rbind(arid1a_clus, arid1b_clus, arid2_clus, multi) 

int_clust <- all %>% 
  group_by(clust, peakname) %>%
  summarise(mean_clust = mean(avg))  %>% 
  filter(mean_clust > 0.10) %>%
  ungroup() %>%
  arrange(desc(mean_clust) )

int_clust$groups <- paste(int_clust$peakname, int_clust$clust, sep = '_') 
dummy <- expand.grid(unique(int_clust$groups), unique(int_clust$groups) )
dummy <- dummy %>% 
  separate(Var1, into = c('startPeak', 'startClust'), sep = '_') %>%
  separate(Var2, into = c('endPeak', 'endClust'), sep = '_') 

# calculate links between clusters
dummy$links <- apply(dummy, 1, createNetworkData, df = all) 
dummy$start <- with(dummy, paste(startPeak, startClust, sep = '_' ) )
dummy$end <- with(dummy, paste(endPeak, endClust, sep = '_') )
# conver to a matrix for network analysis
dx <- dummy %>% 
  select(start, end, links) %>% 
  spread(start, links) 
rownames(dx) <- dx$end
ddx <-dx %>% 
  select(-end) 
diag(ddx) <- 0
# plot the network  aka. fig 5b. 
pal <- c('#fb6a4a','#6baed6', '#74c476', '#756bb1')
cnum <- unlist(sapply(rownames(ddx), function(x) unlist(strsplit(x, split = '_') )[2] ) )
cnames <- unlist(sapply(rownames(ddx), function(x) unlist(strsplit(x, split = '_') ) [1] ) )
clrs <- pal[as.factor(cnames)]
pdf('figures/Fig_5B_network.pdf') 
qgraph(ddx, layout = 'spring', labels =cnum, groups = cnames, color = pal, edge.color = 'grey10', legend.cex = 1.25, GLratio = 1.5)  
dev.off()

# Supplemental Figure 6-2
# ------------------------------------------------------------------------------
mybl2 <- 'Mybl2sc81192V0422111Aln_s_coverage.h5'
mbd4 <- 'Mbd4sc271530V0422111Aln_s_coverage.h5'
tead4 <- 'Tead4sc101184V0422111Aln_s_coverage.h5'
nfic <- 'Nficsc81335V0422111Aln_s_coverage.h5'
p <- plot_arids_by_type(mybl2, 'arid1b_p_', arid1b)
legend <- g_legend(p) 

mybl2_p <- plot_arids_by_type(mybl2, 'arid1b_p_', arid1b) + theme(legend.position = 'none')
mbd4_p <- plot_arids_by_type(mbd4, 'arid1b_p_', arid1b) + theme(legend.position = 'none') + ylab('') 
tead4_p <- plot_arids_by_type(tead4, 'arid1b_p_', arid1b) + theme(legend.position = 'none') + ylab('') 
nfic_p <- plot_arids_by_type(nfic, 'arid1b_p_', arid1b) + theme(legend.position = 'none') + ylab('') 
leg <- plot(legend)

pdf('figures/supplemental/S5_2_metaplots.pdf',  width = 24, 8)
grid.arrange(mybl2_p, mbd4_p, tead4_p, nfic_p, nrow = 1) 
dev.off()

# Show clustering by group
# Figure 6C
# -----------------------------------------------------------------------------=
hmcol <- colorRampPalette(brewer.pal(5, 'RdBu'))(20)
dist.pear <- function(x) as.dist(1-cor(t(x) ) )
hclust.ave <- function(x) hclust(x, method = 'average')
tree <- function(x) cutree(hclust(x, method = 'average'), 3)

plotHM <- function(df) { 
  par(mar = c(0,0,0,0))
  dfm <- df %>% 
   select(EZH2, H3K27ME3, H3K4ME2, H3K27AC, 
          ARID3A, MAX, YY1, BRCA1,  MYBL2, 
          TEAD4, NFIC, MBD4, HDAC2, FOSL2, HNF4A)
  tree_df <- cutree(hclust.ave(dist.pear(as.matrix(dfm))), 3) 
  rowcol <- brewer.pal(3, 'Set1')
  heatmap.2(as.matrix(dfm),
                 distfun = dist.pear, 
                 hclustfun = hclust.ave, 
                 dendrogram = 'none',
                 scale = 'none', 
                 Colv=F,
                 RowSideColors = rowcol[tree_df],
                 margins = c(0,0), 
                 trace = 'none',
                 col = rev(hmcol ), 
                 labRow = F, 
                 labCol = F, 
                 key = F) 
}



# None of these save with the labels
# see the above function for the labels, or check Figure 6C. 
tiff('figures/Fig_5C_arid1a.tiff') 
arid1a_enrich %>% plotHM
dev.off()

tiff('figures/Fig_5C_arid1b.tiff')
arid1b_enrich %>% plotHM
dev.off()

tiff('figures/Fig_5C_arid2.tiff') 
arid2_enrich %>% plotHM
dev.off()

tiff('figures/Fig_5C_multi.tiff')
multi_enrich %>% plotHM
dev.off()

