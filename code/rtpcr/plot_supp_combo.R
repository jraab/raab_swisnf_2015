# Script to combine and plot results from combinatorial KD on second set of siRNAs
library(ggplot2) 
library(dplyr) 
library(broom)
source('code/util/theme_paper.R')

x0182 <- read_csv('data/experimental/rtpcr/val/x0182_combo_results.csv')
x0190 <- read_csv('data/experimental/rtpcr/val/x0190_combo_results.csv')
colnames(x0190) <- c('Primer', 'Sample', 'norm')
x0196 <- read_csv('data/experimental/rtpcr/val/x0196_results_comb.csv') 

comb <- rbind(x0182, x0190, x0196) 

v <- comb %>% 
  mutate(Sample = ifelse(Sample == 'siARID1A_siARID2', 'siARID1A_ARID2', Sample)) %>%
  group_by(Primer, Sample) %>%
  summarise(u  = mean(norm), ci =  (sd(norm)/sqrt(n())) ) %>%
  mutate(Sample = factor(Sample, levels = c('siARID1A', 'siARID2', 'siARID1A_ARID2')))

arids <- v %>%
  filter(!Sample == 'NTG') %>%
  filter(Primer %in% c('ARID1A', 'ARID2'))  %>%
  ggplot(aes(x = Sample, y = u)) + 
  geom_point(size = 4) + 
  geom_errorbar(aes(ymin = u-ci, ymax = u+ci), width = 0.1, size = 2) + 
  geom_hline(y= 0, color = 'red2', linetype = 'dashed') + 
  facet_wrap(~Primer, nrow = 1) + 
  theme_paper() + xlab('') + ylab('Log2 Fold Change') + 
  theme(axis.text = element_text(color = 'grey20'), 
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_line(color= 'grey80', linetype = 'dashed'), 
        panel.grid.major.x = element_blank(), 
        axis.line = element_blank(), 
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(color = 'grey30', fill=NA)) 
ggsave('figures/supplemental/new10_arids.pdf', height = 6, width = 9.6)

tests <- v %>%
  filter(!Sample == 'NTG') %>%
  filter(Primer %in% c('EPHA2', 'RRAS', 'S100A14', 'SPOCK2', 'TNS4'))  %>%
  ggplot(aes(x = Sample, y = u)) + 
  geom_point(size = 4) + 
  geom_errorbar(aes(ymin = u-ci, ymax = u+ci), width = 0.1, size = 2) + 
  geom_hline(y= 0, color = 'red2', linetype = 'dashed') + 
  facet_wrap(~Primer, nrow = 1) + 
  theme_paper() + xlab('') + ylab('Log2 Fold Change') + 
  theme(axis.text = element_text(color = 'grey20'), 
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_line(color= 'grey80', linetype = 'dashed'), 
        panel.grid.major.x = element_blank(), 
        axis.line = element_blank(), 
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(color = 'grey30', fill=NA)) 

ggsave('figures/supplemental/new10_comboval_test.pdf', tests, height = 6, width = 24)

comb %>%
  group_by(Primer) %>%
  filter(!Sample == 'NTG') %>%
  aov(norm ~ Sample, data = .) %>%
  TukeyHSD()

