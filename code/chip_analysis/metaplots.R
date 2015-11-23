# Script to plot metagene profiles over arid peaks smmmits of different classes 
# -----------------------------------------------------------------------------
# 
#   NEED TO COME UP WITH OUT TO USE THE CLUSTERING TO PLOT THE MOST INTERESTING
#   GROUPS USING THIS SCRIPT
#   PROBABLY AS SOME SORT OF HEATMAP Y-ORDERD CONSISTANT 
#   2-3 FACTORS PER GROUP THAT HOPEFULLY VARY SOMEWHAT ACROSS ARIDS
#  See the very bottom if you want to just do a single group. 
# Figure 6 made with this script. - Should clean up and make auto work at some point. 
# -----------------------------------------------------------------------------
library(dplyr)
library(tidyr) 
library(magrittr)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(Repitools)
library(gridExtra)
source('code/theme_paper.R')
#  Set up files and directories
# -----------------------------------------------------------------------------
#PROJ <- '/Users/jraab/proj/swi_snf_final/' # on mac
#PROJ <- '/proj/magnuslb/users/jraab/swi_snf_final/' # on kure
PROJ <- '/home/jraab/proj/swi_snf_final/' # on home
options(stringsAsFactors = FALSE) 
dirad <- 'output/encode_coverages/coverages/'
# also need the arid1a, arid1b, and arid2 alone peaks names - get from annotation
arid1a <- as.list(read.csv('output/macs_peaks/cleaned/arid1a_annotated.csv') %>% 
          filter(svm == 'Single') %>% select(name) )
arid1a_all <- read.csv('output/macs_peaks/cleaned/arid1a_annotated.csv') 

arid1b <- as.list(read.csv('output/macs_peaks/cleaned/arid1b_annotated.csv') %>% 
                    filter(svm == 'Single') %>% select(name) )
arid1b_all <- read.csv('output/macs_peaks/cleaned/arid1b_annotated.csv') 

arid2 <- as.list(read.csv('output/macs_peaks/cleaned/arid2_annotated.csv') %>% 
                   filter(svm == 'Single') %>% select(name) )
arid2_all <- read.csv('output/macs_peaks/cleaned/arid2_annotated.csv') 

multiple <- as.list(read.csv('output/macs_peaks/cleaned/multi_bound_peaks.csv', header =F) %>% 
                      select(V4) ) 
peaks_names <-c(multiple, arid1a, arid1b, arid2) 
# TxnFactors/HistoneMods To look at
#-----------------------------------------------------------------------------
# This is mostly just accounting to organize what I want to do below 
# This will be automated when I've settled on how things will be plotted and
#   organized 
# stubs will control the first half of file name
stubs <- c('multi', 'arid1a', 'arid1b', 'arid2') 
# ends <- c('_p_H3k04me1StdAln_s_coverage.h5', '_p_H3k4me3StdAln_s_coverage.h5', '_p_H3k27acStdAln_s_coverage.h5',
#           '_p_DNase_s_coverage.h5', '_p_Chd2ab68301IggrabAln_s_coverage.h5', '_p_Hdac2sc6296V0416101Aln_s_coverage.h5', 
#           '_p_Arid3anb100279IggrabAln_s_coverage.h5', '_p_Yy1sc281V0416101Aln_s_coverage.h5', '_p_Brca1a300IggrabAln_s_coverage.h5', 
#           '_p_Ezh239875Aln_s_coverage.h5', '_p_Foxa1sc101058V0416101Aln_s_coverage.h5', '_p_Mxi1StdAln_s_coverage.h5', 
#           '_p_MaxIggrabAln_s_coverage.h5', '_p_Hnf4asc8987V0416101Aln_s_coverage.h5', '_p_Hnf4gsc6558V0416101Aln_s_coverage.h5', 
#           '_p_H2azStdAln_s_coverage.h5', '_p_P300sc582IggrabAln_s_coverage.h5', '_p_Fosl2V0416101Aln_s_coverage.h5', 
#           '_p_TbpIggrabAln_s_coverage.h5', '_p_Tead4sc101184V0422111Aln_s_coverage.h5')
# titles <- c('H3K4me1', 'H3K4me3', 'H3K27ac', 'Dnase', 'Chd2', 'Hdac2', 'Arid3a', 'Yy1', 
#             'Brca1', 'Ezh2', 'Foxa1', 'Mxi1', 'Max', 'Hnf4a', 'Hnf4g', 'H2az', 'P300', 
#             'Fosl2', 'Tbp', 'Tead4')
ends <- c('_p_H3k04me1StdAln_s_coverage.h5', '_p_H3k4me3StdAln_s_coverage.h5', '_p_H3k27acStdAln_s_coverage.h5',
           '_p_Arid3anb100279IggrabAln_s_coverage.h5', '_p_Yy1sc281V0416101Aln_s_coverage.h5','_p_Foxa1sc101058V0416101Aln_s_coverage.h5',
          '_p_MaxIggrabAln_s_coverage.h5', '_p_Hnf4gsc6558V0416101Aln_s_coverage.h5', 
          '_p_Fosl2V0416101Aln_s_coverage.h5','_p_Tead4sc101184V0422111Aln_s_coverage.h5', '_p_Mybl2sc81192V0422111Aln_s_coverage.h5', '_p_RxraPcr1xAln_s_coverage.h5')
titles <- c('H3K4me1', 'H3K4me3', 'H3K27ac', 'Arid3a', 'Yy1', 
             'Foxa1', 'Max', 'Hnf4g', 'Fosl2', 'Tead4', 'Mybl2', 'Rxra')




# Libraries
# -----------------------------------------------------------------------------
library(dplyr)
library(ggplot2) 
library(rhdf5) 

# Functions 
# -----------------------------------------------------------------------------
readh5vals <- function(fn) {
  tmp <- h5read(fn, 'coverages', bit64conversion='double')
  vals <- t(tmp$block0_values) 
  row.names(vals) <- tmp$axis1
  return(vals) 
}

readh5vals_fortype <- function(fn) {
  tmp <- h5read(fn, 'coverages', bit64conversion='double')
  vals <- data.frame(t(tmp$block0_values) )
  vals$rowname <- tmp$axis1
  print(dim(vals) )
  return(vals) 
}

filter_arid_coverage <- function(arid_coverage_file, peak_names) { 
  # will return just the coverage for the correct list of peaks
  vals <- readh5vals(arid_coverage_file) 
  return(vals[row.names(vals) %in% peak_names]) 
}

collapsedf <- function(list_dfs) { 
  odf <- data.frame()
  for (i in 1:length(list_dfs) ) { 
    odf <- rbind(odf, list_dfs[[i]])
  }
  return(odf)
}

melt_df_list <- function(df_list) { 
  for (i in 1:length(df_list)){
    df <- df_list[[i]]
    df$Peakname <- row.names(df) 
    mf <- gather(df, 'Position', 'Value', -Peakname) 
  }
  return(output_df) 
}

ci <- function(vals) { 
  return( qnorm(0.975) * sd(vals, na.rm =T) /sqrt(length(vals) ) )
}



plot_four_groups <- function(end, stubs, peaknames) { 
  cols <- c('#ef3b2c', '#67a9cf', '#7fbf7b',  '#737373') 
  vals <- sapply(stubs, function(x) readh5vals(paste0(dirad, x, end)))
  #vals needs filtered
  
  o <- lapply(seq_along(vals),
             function(x, v) { as.data.frame(v[[x]]) %>% mutate(PeakNames = row.names(.) ) %>% 
                                filter(PeakNames %in% peaknames[[x]]) %>%
                                gather('Position', 'Value', -PeakNames) %>% 
                                mutate(group = as.character(rep(names(v[x])) )  ) %>% 
                                group_by(Position, group) %>% 
                                summarise(avg = mean(Value, na.rm =T), ci = ci(Value) ) %>% 
                                mutate (ymin = avg - ci, ymax = avg +ci )}, v = vals) 
  
  output_df <- collapsedf(o) 
  output_df$Position <- rep(seq(-2500, 2490, by = 10), 4) 
  output_df$group <- sapply(output_df$group, simpleCap) 
  g <- ggplot(output_df, aes(x = Position, y = avg, color = group, fill = group) )
  g <- g + geom_line(show_guide = F, size = 1.5) +  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.25, color = NA ) 
  g <- g +theme_minimal() + xlab('') + ylab('IP/Input') 
  g <- g + theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank(), 
                 axis.text = element_text(size = 20), axis.title = element_text(size = 24), 
                 legend.text = element_text(size = 24))  
  g <- g + scale_fill_manual(values = cols) 
  g <- g + scale_color_manual(values = cols) 
  return(g) 
}

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}

plot_arids <- function(end, stubs) { 
    # plot metagene plots for all aird peaks for a given hdf5 files
  cols <- c('#ef3b2c', '#67a9cf', '#7fbf7b',  '#737373') 
  vals <- sapply(stubs, function(x) readh5vals(paste0(dirad, x, end) ) ) 
  o <- lapply(seq_along(vals), 
              function(x, v) { as.data.frame(v[[x]]) %>%  
                              mutate(PeakNames = row.names(.) ) %>% 
                              gather('Position', 'Value', -PeakNames) %>% 
                              mutate(group = as.character(rep(names(v[x] ) ) ) ) %>% 
                              group_by(Position, group) %>% 
                              summarise(avg = mean(Value, na.rm = T), ci = ci(Value) ) %>% 
                              mutate (ymin = avg - ci, ymax = avg + ci)
              }, v = vals )
  output_df <- collapsedf(o) 
  output_df$Position <- rep(seq(-2500, 2490, by = 10), length(stubs)) # only 3 groups this time
  output_df$group <- sapply(output_df$group, simpleCap) 
  g <- ggplot(output_df, aes(x = Position, y = avg, color = group, fill = group) )
  g <- g + geom_line(show_guide = F, size = 1.5) +  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.25, color = NA ) 
  g <- g +theme_minimal() + xlab('') + ylab('IP/Input') 
  g <- g + theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank(), 
                 axis.text = element_text(size = 20), axis.title = element_text(size = 24), 
                 legend.text = element_text(size = 24))  
  g <- g + scale_fill_manual(values = cols) 
  g <- g + scale_color_manual(values = cols) 
  return(g) 
}


plot_arids_by_type <- function(end, stubs, peak_info) { 
  # same as above but plot the arids distinguised by genomic feature

  vals <- lapply(stubs, function(x) readh5vals_fortype(paste0(dirad, x, end) ) ) 
  #vals <- lapply(vals, function(x) as.data.frame(x) )
  #vals <- lapply(vals, function(x) x$rowname <- row.names(x) )
  o <- lapply(seq_along(vals), 
              function(x, v) { v[[x]] %>% tbl_df() %>% 
                  full_join(peak_info[[x]], by = c('rowname'='name')) %>% 
                  select(-nearestTSS, -distance, -svm, -state, -start, -end, -seqnames) %>% 
                  gather('Position', 'Value', -rowname, -type) %>%
                  mutate(group = rep(names(peak_info[x] ) ) )  %>% 
                  filter(!type == 'Exon') %>%
                  group_by(Position, group, type) %>% 
                  summarise(avg = mean(Value, na.rm = T), ci = ci(Value) ) %>% 
                  mutate (ymin = avg - ci, ymax = avg + ci) 
              }, v = vals )
  output_df <- collapsedf(o) 
  
  output_df$Position <- rep(seq(-2500, 2490, by = 10), each = 3, times = 3) # only 3 groups this time
  g <- ggplot(output_df, aes(x = Position, y = avg, color = type , fill = type) )
  g <- g + geom_line(show_guide = F, size = 1.5) +  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.25, color = NA ) 
  g <- g +theme_minimal() + xlab('') + ylab('IP/Input') +theme_paper()
  g <- g + theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank(), 
                 axis.text = element_text(size = 20), axis.title = element_text(size = 24), 
                 legend.text = element_text(size = 24))  
  g <- g + facet_wrap(~group,nrow =1) 

  return(g) 
}

getfn <- function(path, rt_patt) { 
  # return list of files that match the pattern
  patt <- paste0('.*_p_',rt_patt,'.*_s_coverage.h5')
  files <- list.files(path = path, patt = patt, full.names = T, ignore.case = T )
  files <- files[!grepl(pattern = 'snf5*', files)] 
  return(files) 
}

getFilenames <- function(arid, endings) { 
  # return list of files that match the pattern
  path <- 'output/encode_coverages/coverages/'
  file_patt <- lapply(endings, function(x) paste0(arid, x) ) 
  files <- unlist(lapply(file_patt, function(x) list.files(path, x, full.names =T, ignore.case =T) ))
  print (files) 
  return(files) 
}
  
# Body
# -----------------------------------------------------------------------------
# need read in each of the 4 files and filter out the arid peaks that are not alone

# identify 8 mods/txn factors worth looking at and create a panel for each of the four peak groups
arid1a_val <- lapply(getFilenames('arid1a', ends), function(x) readh5vals(x) )
arid1a_val <- lapply(arid1a_val, function(x) x[row.names(x) %in% arid1a[[1]], ] )

arid1b_val <-lapply(getFilenames('arid1b',ends), function(x) readh5vals(x) )
arid1b_val <- lapply(arid1b_val, function(x) x[row.names(x) %in% arid1b[[1]], ])

arid2_val <-lapply(getFilenames('arid2', ends), function(x) readh5vals(x) )
arid2_val <- lapply(arid2_val, function(x) x[row.names(x) %in% arid2[[1]], ])

multi_val <- lapply(getFilenames('multi', ends), function(x) readh5vals(x) )

l <- lapply(ends, function(x, y) plot_four_groups(x, y, c(multiple, arid1a, arid1b, arid2)), y = stubs )

g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

legend <- g_legend(l[[1]])
grid.newpage()
grid.draw(legend)

ll <- sapply(l, function(x) x + theme(legend.position = 'none') ) 

l[[1]] + theme(legend.position = 'none', plot.title = element_text(size = 28))  + ggtitle(titles[1] ) 
l[[2]] + theme(legend.position = 'none', plot.title = element_text(size = 28)) + ggtitle(titles[2]) 
l[[3]] + theme(legend.position = 'none', plot.title = element_text(size =28)) + ggtitle(titles[3]) 
l[[4]] + theme(legend.position = 'none', plot.title = element_text(size = 28)) + ggtitle(titles[4]) 
l[[5]] + theme(legend.position = 'none', plot.title = element_text(size = 28)) + ggtitle(titles[5]) 
l[[6]] + theme(legend.position = 'none', plot.title = element_text(size = 28)) + ggtitle(titles[6]) 
l[[7]] + theme(legend.position = 'none', plot.title = element_text(size = 28)) + ggtitle(titles[7]) 
l[[8]] + theme(legend.position = 'none', plot.title = element_text(size = 28)) + ggtitle(titles[8]) 
l[[9]] + theme(legend.position = 'none', plot.title = element_text(size = 28)) + ggtitle(titles[9]) 
l[[10]] + theme(legend.position = 'none', plot.title = element_text(size = 28)) + ggtitle(titles[10]) 
l[[11]] + theme(legend.position = 'none', plot.title = element_text(size = 28)) + ggtitle(titles[11]) 
l[[12]] + theme(legend.position = 'none', plot.title = element_text(size = 28)) + ggtitle(titles[12])

# H3K4me1
plot_arids(ends[1], stubs[2:4]) + theme(legend.position = 'none') 
# H3K4me3
plot_arids(ends[2], stubs[2:4])
# H3K27ac
plot_arids(ends[3], stubs[2:4]) + theme(legend.position = 'none') 

# can easily check a singly factor down here. 
plot_arids('_p_Sp1Pcr1xAln_s_coverage.h5', stubs)
plot_arids('_p_Tcf7l2UcdAln_s_coverage.h5', stubs)
plot_arids('_p_Hdac2sc6296V0416101Aln_s_coverage.h5', stubs)

arid_dfs <- list(tbl_df(arid1a_all), tbl_df(arid1b_all), tbl_df(arid2_all))
names(arid_dfs) <- c('Arid1a', 'Arid1b', 'Arid2') 
plot_arids_by_type('_p_Hdac2sc6296V0416101Aln_s_coverage.h5', stubs[2:4],  arid_dfs) 
plot_arids_by_type('_p_Sp1Pcr1xAln_s_coverage.h5', stubs[2:4], arid_dfs)
plot_arids_by_type('_p_Tcf7l2UcdAln_s_coverage.h5', stubs[2:4], arid_dfs)
plot_arids_by_type('_p_Hdac2sc6296V0416101Aln_s_coverage.h5', stubs[2:4], arid_dfs)
plot_arids_by_type('_p_Hnf4asc8987V0416101Aln_s_coverage.h5', stubs[2:4], arid_dfs) 
plot_arids_by_type('_p_Hnf4asc8987V0416101Aln_s_coverage.h5', stubs[2:4], arid_dfs) 
plot_arids_by_type('_p_Foxa1sc101058V0416101Aln_s_coverage.h5', stubs[2:4], arid_dfs) 
plot_arids_by_type('_p_Foxa1sc6553V0416101Aln_s_coverage.h5', stubs[2:4], arid_dfs) 
plot_arids_by_type('_p_MaxIggrabAln_s_coverage.h5', stubs[2:4], arid_dfs) 
plot_arids_by_type('_P_Mybl2sc81192V0422111Aln_s_coverage.h5', stubs[2:4], arid_dfs) 

