# Chipseq Figures (2) - to look at histone and txn factor enrichment at peaks 
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
# Important Directories
# -----------------------------------------------------------------------------
PROJ <- '/Users/jraab/proj/swi_snf_final/' # on mac
#PROJ <- '/magnuson-lab/jraab/analysis/swi_snf_final/' # on codon
setwd(PROJ)
output_figs <- paste0(PROJ, 'output/plots/chip/')

# Function 
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
  #d <- hclust(mahalanobis(as.matrix(df_hm), colMeans(df_hm), cov(as.matrix(df_hm)) ), method = 'complete')
  cs_cols <- brewer.pal(5, 'Set1')[as.numeric(getdendro(df_hm, 5) ) ]
  heatmap.2(as.matrix(df_hm), col = rev(hmcol), trace = 'none', scale = 'none', ColSideColors = cs_cols)  
}

cleanDF <- function(df) { 
  droppatt <- 'Pol2Forskln*|Srebp1Insln*|Pgc1aForskln*|Mafk_sc*|Hsf1Forskln*|Grp20Forskln*|ErraForskln*|RxlchPcr1x*|RxlchPcr2x*|RxlchV0416101*|RxlchV0422111*'
  df <- df[!grepl(droppatt, x = row.names(df)), !grepl(droppatt, x = colnames(df) )]
  return(df) 
}



cleanNames <- function(name_vector) { 
  name_vector = gsub(pattern = 'Igg.*$|Ucd.*$|Std.*$|Aln.*$', replacement = '', x=name_vector, perl = T)
  name_vector = gsub(pattern = 'Arid3anb100279', replacement='Arid3a', name_vector, perl=T)
  name_vector = gsub(pattern = 'Bhlhe40c', replacement='Bhlhe40_a', name_vector, perl=T)
  name_vector = gsub(pattern = 'Brca1a300', replacement='Brca1', name_vector, perl=T)
  name_vector = gsub(pattern = 'Chd2ab68301', replacement='Chd2', name_vector, perl=T)
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

#this is not implemented
read_in_data <- function(filename, scale =T){ 
  df <- fread(filename)
  if (scale == TRUE){ 
    e <- scale(df$enrichment, scale =T, center =T) )
  else{ 
    e <- df$enrichment
  }
}
# Load Data
# -----------------------------------------------------------------------------
arid1a <- read.table('output/macs_peaks/cleaned/arid1a_annotated.csv', sep =',', header =T) 
arid1b <- read.table('output/macs_peaks/cleaned/arid1b_annotated.csv', sep =',', header =T)  
arid2 <- read.table('output/macs_peaks/cleaned/arid2_annotated.csv', sep = ',', header =T) 
multi <- read.table('output/macs_peaks/cleaned/multi_bound_peaks.csv', sep = ',', header =F) # multi doesn't seem to have a name
colnames(multi) <- c('seqnames', 'start', 'end','name')  
 
# Get information for enrichment
# ------------------------------------------------------------------------------
enrichment_dir <- 'output/encode_coverages/enrichments/'
a <- read.table(paste0(enrichment_dir, 'multi_p_Znf274UcdAln_s_enrichment.csv'), header =T, sep= ',')

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

# fix the names of the antibodies
# several have other annoying info these have to be dealt with one at a time
###########
df$signame = gsub(pattern = 'Igg.*$|Ucd.*$|Std.*$|Aln.*$', replacement = '', x=df$signame, perl = T)
df$signame = gsub(pattern = 'Arid3anb100279', replacement='Arid3a', df$signame, perl=T)
df$signame = gsub(pattern = 'Bhlhe40c', replacement='Bhlhe40_a', df$signame, perl=T)
df$signame = gsub(pattern = 'Brca1a300', replacement='Brca1', df$signame, perl=T)
df$signame = gsub(pattern = 'Chd2ab68301', replacement='Chd2', df$signame, perl=T)
df$signame = gsub(pattern = '^Ctcf$', replacement = 'Ctcf_a', df$signame, perl = T) 
df$signame = gsub(pattern = 'Corestsc30189', replacement='Corest', df$signame, perl=T)
df$signame = gsub(pattern = 'Ezh239875', replacement='Ezh2', df$signame, perl=T)
df$signame = gsub(pattern = 'Maffm8194', replacement='Maff', df$signame, perl=T)
df$signame = gsub(pattern = 'Mafkab50322', replacement='Mafk_a', df$signame, perl=T)
df$signame = gsub(pattern = 'Mafksc477', replacement='Mafk_b', df$signame, perl=T)
df$signame = gsub(pattern = '^Max$', replacement = 'Max_a', df$signame, perl = T) 
df$signame = gsub(pattern = 'Mazab85725', replacement='Maz', df$signame, perl=T)
df$signame = gsub(pattern = 'P300sc582', replacement='P300_a', df$signame, perl=T)
df$signame = gsub(pattern = '^Rad21$', replacement = 'Rad21_a', df$signame, perl=T)
df$signame = gsub(pattern = 'Rfx5200401194', replacement='Rfx', df$signame, perl=T)
df$signame = gsub(pattern = 'Smc3.*$', replacement = 'Smc3', df$signame, perl=T)
df$signame = gsub(pattern = 'Zeb1.*$', replacement = 'Zeb1', df$signame, perl=T) 
df$signame = gsub(pattern = 'Pcr1x', replacement = '', df$signame, perl=T) 
df$signame = gsub(pattern = 'Pcr2x', replacement = '', df$signame, perl=T) 
df$signame = gsub(pattern = 'Yy1.*$', replacement = 'Yy1', df$signame, perl=T) 
df$signame = gsub(pattern = 'Atf3.*$', replacement = 'Atf3', df$signame, perl=T) 
df$signame = gsub(pattern = 'Bhlhe40V.*$', replacement = 'Bhlhe40_b', df$signame, perl=T)
df$signame = gsub(pattern = 'Cebpd.*$', replacement = 'Cebpd', df$signame, perl=T) 
df$signame = gsub(pattern = 'Creb1.*$', replacement = 'Creb1', df$signame, perl = T) 
df$signame = gsub(pattern = 'Ctcfsc5916V.*$', replacement = 'Ctcf_b', df$signame, perl = T) 
df$signame = gsub(pattern = 'Elf1.*$', replacement = 'Elf1', df$signame, perl=T) 
df$signame = gsub(pattern = 'Fosl2.*$', replacement = 'Fosl2', df$signame, perl = T) 
df$signame = gsub(pattern = 'Foxa1sc101.*$', replacement ='Foxa1_a', df$signame, perl = T) 
df$signame = gsub(pattern = 'Foxa1sc655.*$', replacement = 'Foxa1_b', df$signame, perl = T) 
df$signame = gsub(pattern = 'Hdac2.*$', replacement = 'Hdac2', df$signame, perl =T) 
df$signame = gsub(pattern = 'Hey1.*$', replacement = 'Hey1', df$signame, perl=T) 
df$signame = gsub(pattern = 'Hnf4a.*$', replacement = 'Hnf4a', df$signame, perl = T) 
df$signame = gsub(pattern = 'Hnf4g.*$', replacement = 'Hnf4g', df$signame, perl = T) 
df$signame = gsub(pattern = 'MaxV.*$', replacement = 'Max_b', df$signame, perl = T) 
df$signame = gsub(pattern = 'Mbd4.*$', replacement = 'Mbd4', df$signame, perl=T) 
df$signame = gsub(pattern = 'Mybl2.*$', replacement = 'Mybl2', df$signame, perl = T) 
df$signame = gsub(pattern = 'Nfic.*$', replacement = 'Nfic', df$signame, perl = T) 
df$signame = gsub(pattern = 'Nr2f2.*$', replacement = 'Nr2f2', df$signame, perl = T) 
df$signame = gsub(pattern = 'NrsfV.*$', replacement = 'Nrsf_h', df$signame, perl = T) 
df$signame = gsub(pattern = 'P300V.*$', replacement = 'P300_b', df$signame, perl = T) 
df$signame = gsub(pattern = 'Rad21V.*$', replacement = 'Rad21_b', df$signame, perl = T) 
df$signame = gsub(pattern = 'Sin3a.*$', replacement = 'Sin3a', df$signame, perl = T) 
df$signame = gsub(pattern = 'Sp2.*$', replacement = 'Sp2', df$signame, perl = T) 
df$signame = gsub(pattern = 'Srf.*$', replacement = 'Srf1', df$signame, perl = T) 
df$signame = gsub(pattern = 'Tead4.*$', replacement = 'Tead4', df$signame, perl =T) 
################


#drop antibodies that had treatedments or the one duplicate (mafk_sc) which both antibodies look similiar in my initial look
dropname <- c('Pol2Forskln', 'Srebp1Insln', 'Pgc1aForskln', 'Mafk_sc', 'Hsf1Forskln', 'Grp20Forskln', 'ErraForskln', 'Rxlch', 'RxlchV0416101', 'RxlchV0422111')
df <- df[!df$signame %in% dropname, ]

#create data frames for each arid and remove peaks that would map to multi
arid1a_df <- df[df$peakname == 'arid1a', ]
arid1a_peaks <- read.csv('output/macs_peaks/cleaned/arid1a_annotated.csv')
arid1a_df <- arid1a_df[arid1a_df$name %in% arid1a_peaks[arid1a_peaks$svm == 'Single', 'name'], ]


arid1b_df = df[df$peakname == 'arid1b', ]
arid1b_peaks = read.csv('output/macs_peaks/cleaned/arid1b_annotated.csv')
arid1b_df <- arid1b_df[arid1b_df$name %in% arid1b_peaks[arid1b_peaks$svm == 'Single', 'name'], ]

arid2_df  = df[df$peakname == 'arid2', ]
arid2_peaks = read.csv('output/macs_peaks/cleaned/arid2_annotated.csv') 
arid2_df <- arid2_df[arid2_df$name %in% arid2_peaks[arid2_peaks$svm == 'Single', 'name'], ]

multi_df <- df[df$peakname == 'multi', ]

#combined data frames
full = rbind(arid1a_df, arid1b_df, arid2_df, multi_df)

full_summary <- full %>% 
  group_by(peakname, signame) %>%
  summarise(avg=mean(enrichment) )

# Pull out just the histones
histones = c('H4k20me1', 'H3k9ac', 'H3k79me2', 'H3k4me2', 'H3k4me3', 'H3k36me3', 
             'H3k27me3', 'H3k27ac', 'H3k09me3', 'H3k04me1', 'H2az')


#separate data by histone vs not
histone_df  <- full_summary[full_summary$signame %in% histones, ]
txn_df <- full_summary[!full_summary$signame %in% histones, ] 

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}

histone_df$peakname <- sapply(histone_df$peakname, simpleCap)
txn_df$peakname <- sapply(txn_df$peakname, simpleCap) 

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


histone_df_plot <- makeheat(histone_df)
ggsave(paste0(output_figs, 'histone_enrichment.pdf'), histone_df_plot, width = 4, height = 5)
txn_df_plot <- makeheat(txn_df)
ggsave(paste0(output_figs, 'txn_enrichment.pdf'), txn_df_plot, width = 4, height = 12) 


# Pairwise r within groups. 
# -----------------------------------------------------------------------------
output_final <- list(NA, NA, NA, NA)
for (k in 1:length(groups)){  
  g <- groups[k]
  print(g) 
  files_interest <- getGroupFiles(g, files) 
  # filter out groups int he dropnames 
  files_interest <- files_interest[!files_interest %in% dropname] 
  signal_names <- sapply(files_interest, getSigname) 
  output <- matrix(rep(0), nrow = length(signal_names), ncol = length(signal_names) )
  
  # Ideally I would do something more cleaver than reading this data in a million times like below
  #enrichment_list <- lapply(files_interest, function(x) read_in_data(paste0(enrichment_dir, x)) )
  
  # now pairwise through the enrihcment list
  for (f1 in 1:length(files_interest)) {
    for (f2 in 1:length(files_interest)){ 
      dat1 = fread(paste0(enrichment_dir,files_interest[f1]) )
      dat2 = fread(paste0(enrichment_dir, files_interest[f2]) )
      r = cor(dat1$enrichment, dat2$enrichment, method = 'spearman')
      print(r)
      output[f1,f2] <- r   
    }
  }
  output  <- data.frame(output) 
  colnames(output) <- signal_names
  rownames(output) <- signal_names
  drops <- c('RxlchPcr1xAln', 'RxlchPcr2xAln', 'RxlchV0416101Aln', 'RxlchV0422111Aln')
  output <- output[!row.names(output) %in% drops , !colnames(output) %in% drops, ]
  output_final[[k]] <- output
}

names(output_final) <- groups

# TODO: save the plots
# TODO: clean these plots up - maybe cut the dendrogram 

plot_hierarchy(output_final$arid1a)

plot_hierarchy(output_final$arid1b) 

plot_hierarchy(output_final$arid2) 

plot_hierarchy(output_final$multi)



getdendro <- function(df, ht) {
  df <- cleanDF(df) 
  row.names(df) <- cleanNames(row.names(df))
  colnames(df) <- cleanNames(colnames(df) )
  dendro <- hclust(dist(as.matrix(df), 'manhattan') )
  vals = cutree(dendro, k=ht)
  ord <- dendro$order
  #plot(dendro, labels = vals, hang = -1, col = vals)
  return(list(ord, dendro, vals) )
}

getClusters <- function(df, clust){ 
  hc <- hclust(dist(as.matrix(df), 'manhattan') )
  df$vals <- cutree(hc, k=clust)
  return(df)
}

#put this into the make heat function - order by dendrogram should be optional
# take the full_summary df and plot a single arid ordered by some clustering function 
plot_hm_single <- function(df, arid, all_cors) {
  grid.newpage()
  # all_cors is the list of output_final from the pairwise correlation 
  # in order of groups below - allows getting cluster values 
  groups <- c('arid1a', 'arid1b', 'arid2', 'multi') # to pull data from the output_final list
  sub <- data.frame(df[df$peakname == arid, ])
  #hc <- long_clusters(sub)
  
   
  #dend <- dendro_data(hc)
  mat <- as.matrix(all_cors[[which(arid == groups)]]) 
  dend <- getdendro(mat, 5) 
  clusters <- dend[[3]]
  ords <- sub[dend[[1]],]$signame
  sub$signame <- factor(sub$signame, levels= ords)
  a <- ggplot(sub, aes(x=peakname, y = signame, fill = avg)) + geom_tile(color = 'grey10')
  a <- a + scale_fill_gradient2(low = 'white', high = 'red')
  a <- a +theme(legend.position = 'none', legend.direction = 'horizontal', 
             axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, size = 18, color = 'grey20'), 
             axis.text.y = element_text(size = 14, color = 'grey20'),
             panel.margin = unit(c(1, 0, 1, 0.5), 'lines')) + xlab('') + ylab('')
  a <- a + coord_fixed()  
  b <- ggdendrogram(dend[[2]], rotate = T) + theme(axis.text = element_blank(), 
                                             panel.margin = unit(c(1, 0.5, 1, 0), 'lines' ) )
  c <- grid.arrange(a, b, ncol =2, clip = T, newpage = T) 
  return(c) 
}


newclust <- hclust(dist(output_final$arid1a))$order


long_clusters <- function(df) { 
  sp <- data.frame(spread(df, peakname, avg) )
  row.names(sp) <- sp$signame
  hc <- hclust(dist(as.matrix(sp)[,-1], method = 'mink' ), method = 'complete' )
  return(hc) 
  }
