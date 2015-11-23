library(dplyr) 
library(magrittr) 
library(ggplot2) 
library(rqpcr) # if you need this use - it is a few functions to do 2^dd(ct) on qpcr data)
# library(devtools); install_github('jraab/rqpcr')
source('code/util/theme_paper.R') 
data_raw <- read.csv('data/experimental/rtpcr/raw_combo.csv')
data_raw$Sample <- factor(data_raw$Sample, levels = c('NTG', 'siArid1a', 'siArid1b', 'siArid2', 'siArid1a_Arid2', 'siArid1b_Arid2' ) )
proc_data <- data_raw %>% process_data() %>% ungroup()
  
arids <- c('ARID1A', 'ARID1B', 'ARID2') 
complete <- c('EPHA2', 'SPOCK2', 'TNS4', 'RRAS', 'S100A14') 
cntrl <- 'GUSB'

# Functions 
return_plot <- function(df, genelist) { 
  return(
    df %>% 
    filter(Primer %in% genelist, !Sample %in% 'NTG') %>% 
    ggplot(aes(x= Sample, y = avg) ) + 
      geom_point(size = 2) + 
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.1 )  + 
      geom_hline(y = 0, color = 'red2', linetype = 'dashed')  + 
      facet_wrap(~Primer, nrow = 1 )  + theme_figure() +xlab('') + ylab('Log2 siArid/siNTG') + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = 'grey10'),
            panel.border = element_rect(color = 'grey30', fill =NA), 
            panel.grid.major = element_line(color = 'grey90'), 
            panel.grid.minor = element_blank(), 
            panel.grid.major.x = element_blank(), 
            strip.background = element_rect(color = 'grey30', fill = NA) ) 
  ) 
}

# Process Raw Data
summary <- meta_control(proc_data, cntrl) 
deltas <- summary %>% 
  group_by(Sample, Group) %>% 
  do(delta_ct_rt(., 'meta') ) 

overneg <- deltas %>% 
  group_by(Primer, Group) %>% 
  do(rel_to_sample(., col = 'Sample', sample ='NTG' , log2=T) )

# Set up averages and 95% ci
final <- overneg %>% 
  group_by(Primer, Sample) %>% 
  summarise(avg = mean(norm, na.rm =T),
            ci = 1.96* sd(norm, na.rm =T)/length(norm) ) %>% 
  mutate(ymin = avg-ci, ymax = avg+ci) 
 

# Plot data subseted by primers
arids_plot <- return_plot(final, arids) 
complete_plot <- return_plot(final, complete) 

arids_plot
complete_plot
ggsave('output/plots/arid_rtpcr_plot.png', arids_plot, width = 4.5, height = 3.5, units = 'in')
ggsave('output/plots/compete_rtpcr_plot.png', complete_plot, width = 7.5, height = 3.5, units = 'in') 

