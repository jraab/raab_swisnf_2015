# Script to create Figure 5
library(dplyr)
library(ggplot2)
library(tidyr)
library(broom) 
library(magrittr)
library(VennDiagram)
library(gridExtra) 
source('code/function_defs/chipseq_fig_fun.R')
source('code/util/theme_paper.R')
output_figs <- 'figures/'
chip_peaks <- 'output/macs_peaks/cleaned/'
# Requires files in 'output/dt_matfiles/aridpeaks_* to run. 

# files for comparing Arids
arid1a <- read.csv('output/macs_peaks/cleaned/arid1a_annotated.csv')
arid1b <- read.csv('output/macs_peaks/cleaned/arid1b_annotated.csv')
arid2  <- read.csv('output/macs_peaks/cleaned/arid2_annotated.csv')
all <- rbind(arid1a, arid1b, arid2) 
all %<>% separate(name, into =c('Arid', 'Name'), sep = '_', extra = 'merge') 
all$type <- factor(as.character(all$type), levels = c('Promoter', 'Exon', 'Intron', 'Distal') )
all$Arid <- toupper(all$Arid) 
all %<>% arrange(type) 

# objects for doing Venn Overlaps
arid1a_gr <- bed2gr(paste0(chip_peaks, 'arid1a_peaks_cf.bed'), stranded = F )
arid1b_gr <- bed2gr(paste0(chip_peaks, 'arid1b_peaks_cf.bed'), stranded = F )
arid2_gr  <- bed2gr(paste0(chip_peaks, 'arid2_peaks_cf.bed' ), stranded = F )

# Venn Diagram (Figure 5A) 
# ------------------------------------------------------------------------------
venn <- c('#e34a33','#2b8cbe', '#2ca25f', '#756bb1', '#636363')
pdf('figures/fig3A_venn3way.pdf')
draw3way(arid1a_gr, arid1b_gr, arid2_gr, c('ARID1A', 'ARID1B', 'ARID2'), 
         c(venn[1], venn[2], venn[3]) )
dev.off()

# Arid signal at Arid peaks (Figure 5C, quant of supplement)
# ------------------------------------------------------------------------------
cols <- c('#ef3b2c', '#67a9cf', '#7fbf7b',  '#737373') 
# Function to read in output from deeptools matrix and split 
# into two groups based on the bed file use for deep tools 
# return a data frame with the rowmean of the signal
split_vals <- function(filename) { 
  val_df <- read.table(filename, skip = 1)
  val_df <- val_df[,8:ncol(val_df)]
  middlebits <- ncol(val_df)/2
  middlebits <- (middlebits-50):(middlebits+50) 
  f <- readLines(filename) 
  split_rows <- which(grepl(pattern = '#', f) )
  arid1a_alone <- val_df[1:(split_rows[1]-1), middlebits]
  print(nrow(arid1a_alone) )
  arid1b_alone <- val_df[(split_rows[1]+1):(split_rows[2]-1), middlebits ]
  print(nrow(arid1b_alone) ) 
  arid2_alone  <- val_df[(split_rows[2]+1):(split_rows[3]-1), middlebits]
  print(nrow(arid2_alone) ) 
  multi  <- val_df[(split_rows[3]+1):(nrow(val_df)-1), middlebits]
  print(nrow(multi) )
  # check that I Only have two groups
  if (length(split_rows) != 4) { 
    print('ERROR: Should only be two groups to split') 
  }
  df <- data.frame(Vals = c(rowMeans(arid1a_alone), rowMeans(arid1b_alone), 
                            rowMeans(arid2_alone), rowMeans(multi) ), 
                   Group = c(rep('ARID1A Alone', nrow(arid1a_alone)), 
                             rep('ARID1B Alone', nrow(arid1b_alone)), 
                             rep('ARID2 Alone', nrow(arid2_alone)), 
                             rep('Multiple', nrow(multi)) ) )
  return(df) 
  
}

arid1a_signal <- 'output/dt_matfiles/aridpeaks_arid1a_mat'
arid1b_signal <- 'output/dt_matfiles/aridpeaks_arid1b_mat'
arid2_signal  <- 'output/dt_matfiles/aridpeaks_arid2_mat'
snf5_signal <- 'output/dt_matfiles/aridpeaks_snf5_mat'

arid1a_df <- split_vals(arid1a_signal) %>% mutate(name = rep('ARID1A Signal') )
arid1b_df <- split_vals(arid1b_signal) %>% mutate(name=rep('ARID1B Signal') )
arid2_df  <- split_vals(arid2_signal) %>% mutate(name = rep('ARID2 Signal') ) 
snf5_df  <- split_vals(snf5_signal) %>% mutate(name = rep('SNF5 Signal') ) 

all <- rbind(arid1a_df, arid1b_df, arid2_df, snf5_df) 

g <- ggplot(all, aes(x = Group, y = log2(Vals), fill = Group) )  + geom_boxplot(notch = T) 
g <- g + facet_wrap(~name, ncol = 4) + theme_figure()
g <- g + xlab('') + ylab('Log2 ARID Signal') 
g <- g + theme(panel.background = element_rect(color = 'grey60'), 
               panel.grid.major.x = element_blank(), 
               axis.ticks = element_blank(), 
               axis.text.x = element_blank(), 
               legend.direction = 'horizontal', 
               legend.position = 'bottom')
g <- g + scale_fill_manual(values = cols)
g
ggsave('figures/fig4c_aridsignal_aridpeaks_bxplt.pdf', g, height =6 , width = 9) 

# Figure 5D - distribution by type
# ------------------------------------------------------------------------------
pal1 <- rev(c('#ef3b2c', '#fc9272', '#2171b5', '#6baed6'))
all <- rbind(arid1a, arid1b, arid2) 
all %<>% separate(name, into =c('Arid', 'Name'), sep = '_', extra = 'merge') 
all$type <- factor(as.character(all$type), levels = c('Promoter', 'Exon', 'Intron', 'Distal') )
all$Arid <- toupper(all$Arid) 
all %<>% arrange(type) 
summ_by_svm <- all %>% group_by(Arid, svm, type ) %>% 
  summarise(number = n() ) %>% 
  mutate(percent = number/sum(number) * 100 )
h <- ggplot(summ_by_svm, aes(x=svm, y = percent, fill = type)) + geom_bar(stat = 'identity')
h <- h+ geom_bar(stat = 'identity', color = 'grey30', show_guide = F) 
h <- h + facet_wrap(~Arid, ncol = 3) 
h <- h + theme_figure() + theme(axis.text = element_text(color = 'grey20'), 
                                axis.ticks = element_blank() )
h <- h + scale_fill_manual(values = pal1) + xlab('') + ylab('Percent') 
h
ggsave('figures/fig4D_genomicdist_byclass.pdf', h) 

#Test significance of changes
arid1a_svm <- summ_by_svm[summ_by_svm$Arid == 'ARID1A', ]
arid1b_svm <- summ_by_svm[summ_by_svm$Arid == 'ARID1B', ]
arid2_svm <- summ_by_svm[summ_by_svm$Arid == 'ARID2', ]

chisq_test_groups <- function(df) { 
  df %>% ungroup() %>% select(-Arid, -percent) %>% 
    spread(type, number) %>% select(-svm) %>% chisq.test() 
}


chisq_test_groups(arid1a_svm) 
chisq_test_groups(arid1b_svm) 
chisq_test_groups(arid2_svm)

# Size of multi vs alone peaks - S5_2
all_peaks <- rbind(arid1a, arid1b, arid2)  
width_bxplt <- all_peaks %>% 
  separate(name, into = c('Arid', 'peaknum'), remove =F, sep ='_', extra = 'merge' ) %>% 
  mutate(width = end-start, Arid= toupper(Arid)) %>% 
  ggplot(aes(x=Arid, y = width, fill = svm)) + 
  geom_boxplot(notch = TRUE, position = position_dodge(0.9)) + 
  scale_fill_manual(values = c('grey30', 'steelblue'), name = 'Peak Type') + 
  theme_paper() + 
  xlab('') + ylab('Width of Peak (bp)') + 
  theme(legend.title = element_text())
ggsave('figures/supplemental/figS4_2_width_of_peaks.pdf')  
width_bxplt

