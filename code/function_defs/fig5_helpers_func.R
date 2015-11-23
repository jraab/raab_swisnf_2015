# Figure 6 Helper functions 
library(ggdendro)
library(data.table)
library(dplyr)
library(ggdendro)
library(RColorBrewer)
library(gridExtra) 
library(ggplot2) 
library(tidyr) 
library(qgraph) 
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

long_clusters <- function(df) { 
  sp <- data.frame(spread(df, peakname, avg) )
  row.names(sp) <- sp$signame
  hc <- hclust(dist(as.matrix(sp)[,-1], method = 'mink' ), method = 'complete' )
  return(hc) 
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
  ct <- cutree(hc, k =9) 
  return(list(hc, ct))
  
}

makeDendroPlot <- function(heat, corplot,  hc, clusters) { 
  cluster_pal <- brewer.pal(length(unique(clusters)), 'Set1')
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
    pdf(paste0('figures/supplemental/S5_1_', name, '_clusterplot.pdf'), height = 10, width = 12, pointsize = 18 )
    
    ly <- grid.arrange(pl[[1]], pl[[2]], pl[[3]], pl[[4]], ncol =4, widths = c(2, 10, 0.5, 2) ) 
    print(ly) 
    dev.off()
    write.table(df, paste0('output/supplemental_data/fig5_', name, '_clusterassign.csv'), row.names =F, col.names =T, sep =',')
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

#drop antibodies that had treatedments or the one duplicate (mafk_sc) which both antibodies look similiar in my initial look
dropname <- c('Pol2Forskln', 'Srebp1Insln', 'Pgc1aForskln', 'MAFK_sc', 'Hsf1Forskln', 'Grp20Forskln', 'ErraForskln', 'Rxlch', 'RxlchV0416101', 'RxlchV0422111')

# Functions for Figure 6-2
readh5vals_fortype <- function(fn) {
  tmp <- h5read(fn, 'coverages', bit64conversion='double')
  vals <- data.frame(t(tmp$block0_values) )
  vals$rowname <- tmp$axis1
  print(dim(vals) )
  return(vals) 
}

ci <- function(vals) { 
  return( qnorm(0.975) * sd(vals, na.rm =T) /sqrt(length(vals) ) )
}

g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 



plot_arids_by_type <- function(end, arid, peak_info) { 
  # plots by group for a single Arid
  vals <- readh5vals_fortype(paste0('output/encode_coverages/coverages/', arid, end) )
  output_df <- vals %>%
    full_join(peak_info, by = c('rowname'='name')) %>% 
    select(-nearestTSS, -distance, -state, -start, -end, -seqnames) %>% 
    gather('Position', 'Value', -rowname, -type, -svm) %>% 
    group_by(Position, type, svm) %>% 
    summarise(avg = mean(Value, na.rm = T), ci = ci(Value) ) %>% 
    mutate (ymin = avg - ci, ymax = avg + ci)  %>% 
    filter(! type == 'Exon') 
  vals %>%
    full_join(peak_info, by = c('rowname'='name')) %>% 
    select(-nearestTSS, -distance, -state, -start, -end, -seqnames) %>% 
    gather('Position', 'Value', -rowname, -type, -svm) %>% 
    group_by(Position, type, svm) %>% 
    summarise(count = n() ) %>% print()
  
  
  output_df$Position <- rep(seq(-2500, 2490, by = 10), each = 6) # only 3 groups this time
  g <- ggplot(output_df, aes(x = Position/1000, y = avg, color = type , fill = type) )
  g <- g + geom_line(show_guide = F, size = 1.5) +  
    geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.5, color = NA ) +
    theme_minimal() + xlab('Distance from TSS (kb)') + ylab('IP/Input') +theme_paper() +
    theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank(), 
          axis.text = element_text(size = 20), axis.title = element_text(size = 24), 
          legend.text = element_text(size = 24), 
          strip.text = element_blank() )  
  g <- g + facet_wrap(~svm,ncol=1) 
  return(g) 
}
