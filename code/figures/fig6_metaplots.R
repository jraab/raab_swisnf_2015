# Figure 7 
# Metaplots
library(dplyr)
library(ggplot2) 
library(rhdf5) 
library(gridExtra)
library(gplots)
source('code/util/theme_paper.R')
source('code/function_defs/chipseq_fig_fun.R')

# Need several files located in 
# output/encode_coverages/coverages 
# keep track of these, I'll put just the ones used in the paper up somewhere for download

# Load in files 
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
multiple <- as.list(read.csv('output/macs_peaks/cleaned/multi_bound_peaks.csv', 
                             header =F) %>% select(V4) ) 
peaks_names <-c(multiple, arid1a, arid1b, arid2) 

# Controls first half of file names 
stubs <- c('multi', 'arid1a', 'arid1b', 'arid2') 
# Controls last half of files names
ends <- c('_p_H3k04me1StdAln_s_coverage.h5', 
          '_p_H3k4me3StdAln_s_coverage.h5', 
          '_p_H3k27acStdAln_s_coverage.h5',
          '_p_Arid3anb100279IggrabAln_s_coverage.h5', 
          '_p_Yy1sc281V0416101Aln_s_coverage.h5',
          '_p_Foxa1sc101058V0416101Aln_s_coverage.h5',
          '_p_MaxIggrabAln_s_coverage.h5', 
          '_p_Hnf4gsc6558V0416101Aln_s_coverage.h5', 
          '_p_Tead4sc101184V0422111Aln_s_coverage.h5',
          '_p_Mybl2sc81192V0422111Aln_s_coverage.h5')
titles <- c('H3K4me1', 'H3K4me3', 'H3K27ac', 'ARID3A', 'YY1', 
            'FOXA1', 'MAX', 'HNF4G', 'TEAD4', 'MYBL2')

g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

l <- lapply(ends, function(x, y) plot_four_groups(x, y, c(multiple, arid1a, arid1b, arid2)), y = stubs )

me1 <- l[[1]] + theme(legend.position = 'none', plot.title = element_text(size = 28))  + ggtitle(titles[1] ) 
me3<- l[[2]] + theme(legend.position = 'none', plot.title = element_text(size = 28)) + ggtitle(titles[2]) 
k27ac <- l[[3]] + theme(legend.position = 'none', plot.title = element_text(size =28)) + ggtitle(titles[3]) 
arid3a <- l[[4]] + theme(legend.position = 'none', plot.title = element_text(size = 28)) + ggtitle(titles[4]) 
yy1 <- l[[5]] + theme(legend.position = 'none', plot.title = element_text(size = 28)) + ggtitle(titles[5]) 
foxa1 <- l[[6]] + theme(legend.position = 'none', plot.title = element_text(size = 28)) + ggtitle(titles[6]) 
max <- l[[7]] + theme(legend.position = 'none', plot.title = element_text(size = 28)) + ggtitle(titles[7]) 
hnf4g <- l[[8]] + theme(legend.position = 'none', plot.title = element_text(size = 28)) + ggtitle(titles[8]) 
tead4 <- l[[9]] + theme(legend.position = 'none', plot.title = element_text(size = 28)) + ggtitle(titles[9]) 
mybl2 <- l[[10]] + theme(legend.position = 'none', plot.title = element_text(size = 28)) + ggtitle(titles[10]) 

pdf('figures/Fig6_A_metaplots.pdf', height = 6, width = 30)
grid.arrange(me3, k27ac, max, yy1, arid3a, nrow = 1) 
dev.off()
pdf('figures/Fig6_B_metaplots.pdf', height = 6, width = 30) 
grid.arrange(me1, tead4, foxa1, mybl2, hnf4g, nrow = 1) 
dev.off()

# Legend is not displayed above, but can be printed here
legend <- g_legend(l[[1]])
grid.newpage()
grid.draw(legend)
