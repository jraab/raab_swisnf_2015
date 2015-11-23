library(readr) 
library(ggplot2) 
library(dplyr) 
library(tidyr)
source('code/util/theme_paper.R')

x0173 <- read_csv('data/experimental/rtpcr/val/x0173_results.csv')
x0190 <- read_csv('data/experimental/rtpcr/val/x0190_results.csv')

comb <- rbind(x0173, x0190, x0196) 
comb <- comb %>%
  mutate(Sample = ifelse(Sample == 'siARID1A-1', 'ARID1A_A', as.character(Sample))) %>%
  mutate(Sample = ifelse(Sample == 'siARID1A-2', 'ARID1A_B', as.character(Sample))) %>% 
  mutate(Sample = ifelse(Sample == 'siARID1B-1', 'ARID1B_A', as.character(Sample))) %>%
  mutate(Sample = ifelse(Sample == 'siARID1B-2', 'ARID1B_B', as.character(Sample))) %>%
  mutate(Sample = ifelse(Sample == 'siARID2-1', 'ARID2_A', as.character(Sample))) %>%
  mutate(Sample = ifelse(Sample == 'siARID2-2', 'ARID2_B', as.character(Sample))) %>%
  separate(Sample, into = c('ARID', 'siRNA'), extra ='drop') %>%
  mutate(siRNA = ifelse(siRNA == 'A', 1, as.character(siRNA))) %>%
  mutate(siRNA = ifelse(siRNA == 'B', 2, as.character(siRNA)))

all_sum <- comb %>% 
  filter(!ARID == 'NTG') %>%
  group_by(Primer, ARID, siRNA) %>%
  summarise(mean_lfc = mean(norm), se = sd(norm)/sqrt(n()))%>% 
  mutate(ymin = mean_lfc - se, ymax = mean_lfc + se) 

cols <- c('grey20', 'grey60')

arid1a <- all_sum %>%
  filter(ARID == 'ARID1A', !Primer %in% c('GUSB', 'meta')) %>%
  ggplot(aes(x = Primer, y = mean_lfc, fill = factor(siRNA))) + 
  geom_errorbar(aes(ymin = ymin, ymax = ymax), position = position_dodge(0.9), width = 0.1) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_fill_manual(values = cols) + theme_paper() + ylab('Log2 Fold Change') + xlab('') + 
  theme(axis.text = element_text(color = 'grey10'))

arid1b <- all_sum %>%
  filter(ARID == 'ARID1B', !Primer %in% c('GUSB', 'meta')) %>% 
  ggplot(aes(x = Primer, y = mean_lfc, fill = factor(siRNA))) + 
  geom_errorbar(aes(ymin = ymin, ymax = ymax), position = position_dodge(0.9), width = 0.1) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_fill_manual(values = cols)  + theme_paper() + ylab('Log2 Fold Change') + xlab('') + 
  theme(axis.text = element_text(color = 'grey10'))

arid2 <- all_sum %>%
  filter(ARID == 'ARID2', !Primer %in% c('GUSB', 'meta')) %>% 
  ggplot(aes(x = Primer, y = mean_lfc, fill = factor(siRNA))) + 
  geom_errorbar(aes(ymin = ymin, ymax = ymax), position = position_dodge(0.9), width = 0.1) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_fill_manual(values = cols)  + theme_paper() + ylab('Log2 Fold Change') + xlab('') + 
  theme(axis.text = element_text(color = 'grey10'))

ggsave('output/diffexp/plots/arid1a_val.pdf', arid1a) 
ggsave('output/diffexp/plots/arid1b_val.pdf', arid1b) 
ggsave('output/diffexp/plots/arid2_val.pdf', arid2) 
