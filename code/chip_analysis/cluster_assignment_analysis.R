# create heatmaps of histonemods/txn factors ordered by *_cor.tab clustering
# ---------------------------------------------------------------------------.
# Script is meant to generate heatmaps of the correlations between different t
# transcription factors and histone modifications. 
# Specifically it will create the plots used in figure and supplemental 5-1
# These include a network chart of the relationship between clusters of the 
# different ARID bound classes, as well as the enrichment and correlation matrix
# plots that show how these clusters were generated. 
#
# This script requires the output of several others
#    pairwise_r_enrichments.R
#    peakannotations (find acutal script that makes the *_annotated.csv files) 
#    coverageFrame.py #run on codon and giving a .csv file for each arid, signal combo
#    
# ---------------------------------------------------------------------------
#Libraries
library(data.table)
library(dplyr)
library(ggdendro)
library(RColorBrewer)
library(gridExtra) 
library(ggplot2) 
library(tidyr) 
library(qgraph) 
#
#
# Impt Directories 
#-----------------------------------------------------------------------------
output_figs <-  'output/plots/chip/cluster_analysis/'
if (!file.exists(output_figs) ){ 
  dir.create(output_figs)
}

# Function Defs
# ----------------------------------------------------------------------------
getGroupFiles <- function(group, filelist){ 
  filelist_groups <- sapply(filelist, function(x) unlist(strsplit(x, split = '_'))[1] )
  keep <- filelist_groups == group
  print(table(keep) )
  return(filelist[keep]) 
}

getSigname <- function(filename){ 
  return(unlist(strsplit(filename, split = '_'))[3])
}

plot_hierarchy <- function(df) { 
  hmcol <- colorRampPalette(brewer.pal(9, 'RdBu'))(100)
  df_hm <- cleanDF(df) 
  colnames(df_hm) <- cleanNames(colnames(df_hm) ) 
  row.names(df_hm) <- cleanNames(row.names(df_hm) )
  heatmap.2(as.matrix(df_hm), col = rev(hmcol), trace = 'none', scale = 'none') 
}

cleanDF <- function(df) { 
  droppatt <- 'Pol2_aForskln|Pol2Forskln*|Srebp1Insln*|Pgc1aForskln*|Mafk_sc*|Hsf1Forskln*|Grp20Forskln*|ErraForskln*|RxlchPcr1x*|RxlchPcr2x*|RxlchV0416101*|RxlchV0422111*'
  df <- df[!grepl(droppatt, x = row.names(df)), !grepl(droppatt, x = colnames(df) )]
  return(df) 
}



cleanNames <- function(name_vector) { 
  name_vector = gsub(pattern = 'Igg.*$|Ucd.*$|Std.*$|Aln.*$', replacement = '', x=name_vector, perl = T)
  name_vector = gsub(pattern = 'Arid3anb100279', replacement='Arid3a', name_vector, perl=T)
  name_vector = gsub(pattern = 'Bhlhe40c', replacement='Bhlhe40_a', name_vector, perl=T)
  name_vector = gsub(pattern = 'Brca1a300', replacement='Brca1', name_vector, perl=T)
  name_vector = gsub(pattern = '^Ctcf$', replacement = 'Ctcf_a', name_vector, perl = T) 
  name_vector = gsub(pattern = 'Corestsc30189', replacement='Corest', name_vector, perl=T)
  name_vector = gsub(pattern = 'Ezh239875', replacement='Ezh2', name_vector, perl=T)
  name_vector = gsub(pattern = 'Maffm8194', replacement='Maff', name_vector, perl=T)
  name_vector = gsub(pattern = 'Mafkab50322', replacement='Mafk_a', name_vector, perl=T)
  name_vector = gsub(pattern = 'Mafksc477', replacement='Mafk_b', name_vector, perl=T)
  name_vector = gsub(pattern = '^Max$', replacement = 'Max_a', name_vector, perl = T) 
  name_vector = gsub(pattern = 'Mazab85725', replacement='Maz', name_vector, perl=T)
  name_vector = gsub(pattern = 'P300sc582', replacement='P300_a', name_vector, perl=T)
  name_vector = gsub(pattern = '^Rad21$', replacement = 'Rad21_a', name_vector, perl=T)
  name_vector = gsub(pattern = 'Rfx5200401194', replacement='Rfx', name_vector, perl=T)
  name_vector = gsub(pattern = 'Smc3.*$', replacement = 'Smc3', name_vector, perl=T)
  name_vector = gsub(pattern = 'Zeb1.*$', replacement = 'Zeb1', name_vector, perl=T) 
  name_vector = gsub(pattern = 'JundPcr1x*', replacement = 'Jund_b', name_vector, perl =T)
  name_vector = gsub(pattern = 'Jund*', replacement = 'Jund_a', name_vector, perl =T) 
  name_vector = gsub(pattern = 'Pol2Pcr2x*', replacement = 'Pol2_b', name_vector, perl =T )
  name_vector = gsub(pattern = 'Pol2*', replacement = 'Pol2_a', name_vector, perl = T)
  name_vector = gsub(pattern = 'Pcr1x', replacement = '', name_vector, perl=T) 
  name_vector = gsub(pattern = 'Pcr2x', replacement = '', name_vector, perl=T) 
  name_vector = gsub(pattern = 'Yy1.*$', replacement = 'Yy1', name_vector, perl=T) 
  name_vector = gsub(pattern = 'Atf3.*$', replacement = 'Atf3', name_vector, perl=T) 
  name_vector = gsub(pattern = 'Bhlhe40V.*$', replacement = 'Bhlhe40_b', name_vector, perl=T)
  name_vector = gsub(pattern = 'Cebpd.*$', replacement = 'Cebpd', name_vector, perl=T) 
  name_vector = gsub(pattern = 'Creb1.*$', replacement = 'Creb1', name_vector, perl = T) 
  name_vector = gsub(pattern = 'Ctcfsc5916V.*$', replacement = 'Ctcf_b', name_vector, perl = T) 
  name_vector = gsub(pattern = 'Elf1.*$', replacement = 'Elf1', name_vector, perl=T) 
  name_vector = gsub(pattern = 'Fosl2.*$', replacement = 'Fosl2', name_vector, perl = T) 
  name_vector = gsub(pattern = 'Foxa1sc101.*$', replacement ='Foxa1_a', name_vector, perl = T) 
  name_vector = gsub(pattern = 'Foxa1sc655.*$', replacement = 'Foxa1_b', name_vector, perl = T) 
  name_vector = gsub(pattern = 'Hdac2.*$', replacement = 'Hdac2', name_vector, perl =T) 
  name_vector = gsub(pattern = 'Hey1.*$', replacement = 'Hey1', name_vector, perl=T) 
  name_vector = gsub(pattern = 'Hnf4a.*$', replacement = 'Hnf4a', name_vector, perl = T) 
  name_vector = gsub(pattern = 'Hnf4g.*$', replacement = 'Hnf4g', name_vector, perl = T) 
  name_vector = gsub(pattern = 'MaxV.*$', replacement = 'Max_b', name_vector, perl = T) 
  name_vector = gsub(pattern = 'Mbd4.*$', replacement = 'Mbd4', name_vector, perl=T) 
  name_vector = gsub(pattern = 'Mybl2.*$', replacement = 'Mybl2', name_vector, perl = T) 
  name_vector = gsub(pattern = 'Nfic.*$', replacement = 'Nfic', name_vector, perl = T) 
  name_vector = gsub(pattern = 'Nr2f2.*$', replacement = 'Nr2f2', name_vector, perl = T) 
  name_vector = gsub(pattern = 'NrsfV.*$', replacement = 'Nrsf_h', name_vector, perl = T) 
  name_vector = gsub(pattern = 'P300V.*$', replacement = 'P300_b', name_vector, perl = T) 
  name_vector = gsub(pattern = 'Rad21V.*$', replacement = 'Rad21_b', name_vector, perl = T) 
  name_vector = gsub(pattern = 'Sin3a.*$', replacement = 'Sin3a', name_vector, perl = T) 
  name_vector = gsub(pattern = 'Sp2.*$', replacement = 'Sp2', name_vector, perl = T) 
  name_vector = gsub(pattern = 'Srf.*$', replacement = 'Srf1', name_vector, perl = T) 
  name_vector = gsub(pattern = 'Tead4.*$', replacement = 'Tead4', name_vector, perl =T) 
  return(name_vector)
}


simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}

makeheat <- function(df) { 
  col= brewer.pal(9,'Reds')
  theme_clean = function() { 
    theme(text=element_text(family='Helvetica'),
          axis.ticks = element_blank(), 
          axis.text.x=element_text(size=16, angle=90, vjust=0.5), 
          axis.text.y=element_text(size=14), 
          legend.title=element_blank(), 
          panel.background=element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(size=16)) 
  }
  p <- ggplot(df, aes(x=peakname, y=signame, fill=avg) ) + geom_tile(colour='grey30') 
  p <- p + theme_clean() + xlab('') + ylab('') 
  p <- p + scale_fill_gradient2(low='white', high='red')  
  return(p) 
}

getClusterInfo <- function(cor_df) { 
  # take the correlation data frame
  # return a list of dendrogram, class assignemnts
  #    Note dendrogram is a proper hclust object, includes labels and order
  # No control over clustering method, or number of classes 
  # Use k = 5  for assignments
  # dendrogram will be cut to k = 5 \
  # method will be DEFAULT for now, but should change
  d <- dist(as.matrix(cor_df)) 
  hc <- hclust(d, method = 'complete')
  ct <- cutree(hc, k =5) 
  return(list(hc, ct))
  
}

makeDendroPlot <- function(heat, corplot,  hc, clusters) { 
  cluster_pal <- brewer.pal(5, 'Set1')
  heat$signame <- factor(heat$signame, levels = heat$signame[hc$order])
  heat$clusters <- factor(clusters)
  corplot <- corplot %>% add_rownames() 
  corplot$signame <- factor(corplot$rowname,  levels = corplot$rowname[hc$order])
  cpm <- gather(corplot, 'signame', 'value', -rowname) 
  cpm$value <- as.numeric(cpm$value) 
  cpm$signame <- factor(cpm$signame, levels = corplot$rowname[hc$order])
  cpm$rowname <- factor(cpm$rowname, levels = corplot$rowname[hc$order]) 
  
  # heatmap of enrichments
  hp <- ggplot(heat, aes(x=1, y = signame)) + geom_tile(color = 'grey10', aes(fill = avg) ) 
  hp <- hp + scale_fill_gradient2(low = 'white', high = '#31a354')  + xlab('') + ylab('') 
  hp <- hp + theme_minimal()
  hp <- hp + theme(axis.ticks = element_blank(), legend.position = 'none', plot.margin = unit(c(1, 0, 1, 1), 'lines')) 
  hp <- hp + theme(axis.text.x = element_blank(), plot.background = element_blank() ) + scale_y_discrete(expand = c(0,0) )
  
  # cluster assignments
  cl <- ggplot(heat, aes(x=1, y = signame, fill = clusters)) + geom_tile()
  cl <- cl + theme_minimal() + theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), legend.position = 'none')
  cl <- cl + theme(plot.margin = unit(c(1, 0, 1.25, 0), 'lines') ) + scale_fill_manual(values = cluster_pal)
  
  # heatmap of correlations
  cpm <- cpm %>% filter(!signame == 'NA', !rowname == 'NA') 
  cp <- ggplot(cpm, aes(x = signame, y = rowname, fill = value)) 
  cp <- cp + geom_tile() 
  cp <- cp + theme_minimal() 
  cp <- cp + theme(axis.ticks = element_blank() ) 
  cp <- cp + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', limits = c(-1,1) )
  cp <- cp + theme(axis.text = element_blank(), 
                   axis.title = element_blank(), 
                   legend.position = 'none', 
                   plot.margin = unit(c(1,0,1.25, 0), 'lines') )
  
  # dendrogram
  dend <- dendro_data(hc, type = 'rectangle')
  dp <- ggplot() + geom_segment(data = segment(dend), aes(x=x, y=y, xend=xend, yend=yend)) + theme_dendro()
  dp <- dp + theme(legend.position = 'none', axis.text = element_blank(), plot.margin = unit(c(1, 1, 1, 0), 'lines') )
  dp <- dp + coord_flip() + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0) )
  
  #return plots above in that order
  return(list(hp, cp, cl, dp) ) 
}

put_together <- function(cormat, heatvals) { 
  clust_info <- getClusterInfo(cormat)
  p <- makeDendroPlot(heatvals, cormat, clust_info[[1]], clust_info[[2]])
  df <- data.frame(signame = names(clust_info[[2]]), clust = clust_info[[2]])
  df <- merge(heatvals, df, by = 'signame') 
  return(list(p, df) ) 
}

save_info <- function(corr, summ, name, print =FALSE) { 
  ret <- put_together(corr, summ) 
  pl <- ret[[1]]
  df <- ret[[2]]
  if (print) { 
    print(grid.arrange(pl[[1]], pl[[2]], pl[[3]], pl[[4]], ncol = 4, widths = c(2, 10, 0.5, 2) ) ) 
  }
  else{
    pdf(paste0(output_figs, name, '_clusterplot.pdf'), height = 10, width = 12, pointsize = 18 )
    
    ly <- grid.arrange(pl[[1]], pl[[2]], pl[[3]], pl[[4]], ncol =4, widths = c(2, 10, 0.5, 2) ) 
    print(ly) 
    dev.off()
    write.table(df, paste0(output_figs, name, '_clusterassign.csv'), row.names =F, col.names =T, sep =',')
    write.table(corr, paste0(output_figs, name, '_corvalues.csv'), row.names = T, col.names = T, sep = ',') 
    return(NULL) 
  }
}
createNetworkData <- function(row, df) { 
  startSigs <- df %>% 
    filter(peakname == row[1] & clust == row[2]) %>% 
    select(signame) %>% mutate(signame = as.character(signame) ) %>% as.list() 
  endSigs <- df %>% 
    filter(peakname == row[3] & clust == row[4]) %>% 
    select(signame) %>% mutate(signame = as.character(signame) ) %>% as.list() 
  retval <- sum(startSigs$signame %in% endSigs$signame) 
  return(retval) 
  
}

subset_enrich_plot <-function(df) { 
  maxscale <- max(all$avg) #bad to call from global env
  g <-  df %>% ggplot(aes(x = peakname, y =signame, fill = avg)) + geom_tile(color = 'grey30') + 
    theme(axis.text = element_text(size = 10, color = 'grey10'), 
          legend.position = 'none', 
          panel.background = element_blank(), 
          plot.background = element_blank(), 
          axis.text.x = element_text(angle = 90 , hjust = 1, vjust = 0.5), 
          axis.ticks = element_blank(), 
          panel.grid = element_blank()) + 
    scale_fill_gradient2(low = 'white', high = '#31a354', limits = c(0,maxscale)) + coord_fixed() + 
    xlab('') + ylab('')
  return(g) 
}

# End function defs
# ------------------------------------------------------------------------------


# Load Data - Information about alonve vs multi peaks
# --------------------------------------------------------------------------
arid1a <- read.table('output/macs_peaks/cleaned/arid1a_annotated.csv', sep =',', header =T) 
arid1b <- read.table('output/macs_peaks/cleaned/arid1b_annotated.csv', sep =',', header =T)  
arid2 <- read.table('output/macs_peaks/cleaned/arid2_annotated.csv', sep = ',', header =T) 
multi <- read.table('output/macs_peaks/cleaned/multi_bound_peaks.csv', sep = ',', header =F) # multi doesn't seem to have a name
colnames(multi) <- c('seqnames', 'start', 'end','name')  

enrichment_dir <- 'output/encode_coverages/enrichments/'

#get enrichment information
files = list.files(enrichment_dir)
groups = c('arid1a', 'arid1b', 'arid2', 'multi') 
df = data.frame(gene=NULL, enrichment=NULL, peakname=NULL, signame=NULL)
for (g in groups){ 
  for (f in files) { 
    peak = unlist(strsplit(f, '_'))[1]
    if (peak %in% g) { 
      dat = fread(paste0(enrichment_dir,f) )
      sig = unlist(strsplit(f, '_'))[3]
      dat$peakname = rep(peak, nrow(dat))
      dat$signame  = rep(sig,  nrow(dat))
      dat$enrichment= as.numeric(as.character(dat$enrichment))
      dat <- dat %>% mutate(enrichment_scaled = (enrichment-min(enrichment))/max(enrichment) ) # scale enrichment 0-1
      df = rbind(df, dat)
    }
  }
}
setnames(df, 'gene', 'name')
df$signame <- cleanNames(df$signame) 

#drop antibodies that had treatments or the one duplicate (mafk_sc) which both antibodies look similiar in my initial look
dropname <- c('Pol2_aForskln', 'Pol2Forskln', 'Srebp1Insln', 'Pgc1aForskln', 'Mafk_sc', 'Hsf1Forskln', 'Grp20Forskln', 'ErraForskln', 'Rxlch', 'RxlchV0416101', 'RxlchV0422111')
df <- df[!df$signame %in% dropname, ]
df$peakname <- sapply(df$peakname, simpleCap) 
df$signame <- sapply(df$signame, simpleCap)

#create data frames for each arid and remove peaks that would map to multi
arid1a_df <- df[df$peakname == 'Arid1a', ]
arid1a_peaks <- read.csv('output/macs_peaks/cleaned/arid1a_annotated.csv')
arid1a_df <- arid1a_df[arid1a_df$name %in% arid1a_peaks[arid1a_peaks$svm == 'Single', 'name'], ]

arid1b_df = df[df$peakname == 'Arid1b', ]
arid1b_peaks = read.csv('output/macs_peaks/cleaned/arid1b_annotated.csv')
arid1b_df <- arid1b_df[arid1b_df$name %in% arid1b_peaks[arid1b_peaks$svm == 'Single', 'name'], ]

arid2_df  = df[df$peakname == 'Arid2', ]
arid2_peaks = read.csv('output/macs_peaks/cleaned/arid2_annotated.csv') 
arid2_df <- arid2_df[arid2_df$name %in% arid2_peaks[arid2_peaks$svm == 'Single', 'name'], ]

multi_df <- df[df$peakname == 'Multi', ]

#combine data frames
full = rbind(arid1a_df, arid1b_df, arid2_df, multi_df)
full_summary <- full %>% 
  group_by(peakname, signame) %>%
  summarise(avg=mean(enrichment) )

arid1a_sum <- full_summary[full_summary$peakname == 'Arid1a', ]
arid1b_sum <- full_summary[full_summary$peakname == 'Arid1b', ]
arid2_sum <- full_summary[full_summary$peakname == 'Arid2',] 
multi_sum <- full_summary[full_summary$peakname == 'Multi', ]

# Now the arid/multi df have information for a single peak set
# and the full_summary contains info for all peaks

# Next use the correlation information from pairwise_r_enrichments.R to cluster signals based on their pairwise similiarites
arid1a_cor <- read.table('output/arid1a_cor.tab', sep ='\t')
arid1b_cor <- read.table('output/arid1b_cor.tab', sep= '\t')
arid2_cor  <- read.table('output/arid2_cor.tab',  sep='\t')
multi_cor  <- read.table('output/multi_cor.tab', sep = '\t')
# for each of those I need a dendrogram, cluster assignments, labels, and order of assignments 


# This saves the plots for figure 5-1 
save_info(arid1a_cor, arid1a_sum, 'arid1a') 
save_info(arid1b_cor, arid1b_sum, 'arid1b') 
save_info(arid2_cor, arid2_sum, 'arid2') 
save_info(multi_cor, multi_sum, 'multi') 

# Uncomment if you runt his manually and just want to see what the plots would look like
save_info(arid1a_cor, arid1a_sum, 'arid1a', print =TRUE) 
save_info(arid1b_cor, arid1b_sum, 'arid1b', print = TRUE) 
save_info(arid2_cor, arid2_sum, 'arid2', print = T) 
save_info(multi_cor, multi_sum, 'multi', print =T) 

# Create a network plot showing relationship between clusters. 
#Share category
arid1a_clus <- read.csv('output/plots/chip/cluster_analysis/arid1a_clusterassign.csv') 
arid1b_clus <- read.csv('output/plots/chip/cluster_analysis/arid1b_clusterassign.csv') 
arid2_clus <- read.csv('output/plots/chip/cluster_analysis/arid2_clusterassign.csv') 
multi <- read.csv('output/plots/chip/cluster_analysis/multi_clusterassign.csv') 

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
head(dummy)
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
pdf('output/plots/chip/cluster_analysis/network.pdf') 
qgraph(ddx, layout = 'spring', labels =cnum, groups = cnames, color = pal, edge.color = 'grey10', legend.cex = 1.5)  
dev.off()

# Pull out the factors invovled in specific networks and makes plots of enrichment 
# This ends up being figure 5c. 
arid1a_c1 <- all %>% 
  filter(peakname == 'Arid1a', clust == 1) 

arid1b_c1 <- all %>% 
  filter(peakname == 'Arid1b', clust == 1) 

arid2_c1 <- all %>% 
  filter(peakname == 'Arid2', clust == 1) 

multi_c1 <- all %>%
  filter(peakname == 'Multi', clust == 1) 

multi_c2 <- all %>% 
  filter(peakname == 'Multi', clust == 2) 

arid1a_c3 <- all %>% 
  filter(peakname == 'Arid1a', clust == 3) 

arid1b_c3 <- all %>% 
  filter(peakname == 'Arid1b', clust == 3 )

arid2_c3 <- all %>% 
  filter(peakname == 'Arid2', clust == 3 ) 

multi_c4 <- all %>% 
  filter(peakname == 'Multi', clust == 4 ) 

arid1b_c5 <- all %>% 
  filter(peakname == 'Arid1b', clust == 5 )

multi_c2 <- all %>% 
  filter(peakname == 'Multi', clust == 2) 

arid1a_c2 <- all %>%
  filter(peakname == 'Multi', clust == 2) 

arid2_c1 <- all %>% 
  filter(peakname == 'Arid2', clust ==1) 

#Using the all data get the names of the signals in each cluster that may be important.

txn_regulators <- union(arid1a_c1$signame, union(arid2_c1$signame, multi_c1$signame) )

txn_reg_plot <- all %>% 
  filter(signame %in% txn_regulators) %>% 
  subset_enrich_plot()
txn_reg_plot
ggsave('output/plots/chip/cluster_analysis/txn_reg_plot.png', txn_reg_plot, dpi = 300, bg = 'transparent')


arid1b_alone <- arid1b_c5$signame

arid1b_alone_plot <- all %>% 
  filter(signame %in% arid1b_alone) %>% 
  subset_enrich_plot() + theme(axis.text.y = element_text(size = 32), axis.text.x = element_blank() )
arid1b_alone_plot
ggsave('output/plots/chip/cluster_analysis/arid1b_alone_plot.png', arid1b_alone_plot, dpi = 300, bg = 'transparent') 

histones <- union(arid1a_c3$signame, union(arid2_c3$signame, multi_c4$signame) ) 

histone_plot <- all %>%
  filter(signame %in% histones) %>%
  subset_enrich_plot() 
histone_plot <- histone_plot + theme(legend.position = 'right', legend.title = element_blank(), axis.text = element_text(size = 20)  ) 
ggsave('output/plots/chip/cluster_analysis/histone_plot.png', histone_plot, dpi = 300, bg = 'transparent') 

# This group is unsaved for now
mystery <- union(arid1a_c2$signame, multi_c2) 
mystery_plot <- all %>%
  filter(signame %in% mystery) %>% 
  subset_enrich_plot() 
mystery_plot <- mystery_plot + theme(plot.background = element_rect(fill = 'transparent', color = 'NA'), 
                     panel.background = element_rect(fill = 'transparent', color = 'NA') ) 

# Done