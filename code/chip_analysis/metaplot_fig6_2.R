# Supplemental Figure 6-2
# Tead4, Mybl2, Nfic, and Mbd4 at Arid1b peaks by location and type
library(ggplot2) 
library(dplyr) 
library(RColorBrewer)
library(rhdf5)
library(tidyr)
source('code/theme_paper.R')
dirad <- 'output/encode_coverages/coverages/'

pal <- brewer.pal(4, 'Set2') 

arid1b <- as.list(read.csv('output/macs_peaks/cleaned/arid1b_annotated.csv') %>% 
                    filter(svm == 'Single') %>% select(name) )
arid1b_all <- read.csv('output/macs_peaks/cleaned/arid1b_annotated.csv') 

# Functions
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
  vals <- readh5vals_fortype(paste0(dirad, arid, end) )
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
  g <- g + geom_line(show_guide = T, size = 1.5) +  
    geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.25, color = NA,show_guide=F ) +
    theme_minimal() + xlab('Distance from TSS (kb)') + ylab('IP/Input') +theme_paper() +
    guides(color = guide_legend(override.aes= list(size =4 ) ) ) + 
    theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank(), 
          axis.text = element_text(size = 20), axis.title = element_text(size = 24), 
          legend.text = element_text(size = 24))  
  g <- g + facet_wrap(~svm,ncol=1) 
  return(g) 
}

###############################################################################
mybl2 <- 'Mybl2sc81192V0422111Aln_s_coverage.h5'
mbd4 <- 'Mbd4sc271530V0422111Aln_s_coverage.h5'
tead4 <- 'Tead4sc101184V0422111Aln_s_coverage.h5'
nfic <- 'Nficsc81335V0422111Aln_s_coverage.h5'
p <- plot_arids_by_type(mybl2, 'arid1b_p_', arid1b_all)
legend <- g_legend(p) 

plot_arids_by_type(mybl2, 'arid1b_p_', arid1b_all) + theme(legend.position = 'none')
plot_arids_by_type(mbd4, 'arid1b_p_', arid1b_all) + theme(legend.position = 'none')
plot_arids_by_type(tead4, 'arid1b_p_', arid1b_all) + theme(legend.position = 'none')
plot_arids_by_type(nfic, 'arid1b_p_', arid1b_all) + theme(legend.position = 'none')
plot(legend)
