library(ggplot2) 
library(dplyr)
library(tidyr)
library(broom)
library(readr) 
source('code/util/theme_paper.R')
comb <- read_csv('data/experimental/chip_kd_results.csv') 

comb <- comb %>% 
  separate(Line, into = c('Line', 'extra'), extra = 'drop') %>%
  filter(!Primer %in% c('RRAS_2', 'VASP_2', 'TNFSRF12A_2'))%>%
  separate(Primer, into =c('Primer', 'extra2'), extra = 'drop') %>%
  mutate(Primer = ifelse(Primer == 'TNFRS12', 'TNFRSF12A', as.character(Primer))) %>%
  mutate(Primer = ifelse(Primer == 'TNFRS12A', 'TNFRSF12A', as.character(Primer))) %>%
  mutate(Line = ifelse(Line == 'NS', 'shNS', as.character(Line)), 
         Primer = factor(as.character(Primer), 
                         levels = c('INS', 'GLOB', 'EPHA2', 'GPSM3', 
                                    'JUN', 'NKD1', 'PPL', 'RRAS', 'TNFRSF12A',
                                    'VASP')) )

comb_summ <- comb %>% 
  group_by(Line, Primer, Sample) %>%
  summarise(u = mean(norm_percent), se = sd(norm_percent)/sqrt(n())) 


plt <- comb_summ %>%
  ungroup() %>%
  filter(Sample %in% c('ARID1A', 'ARID2', 'H3K27ac'), 
         Primer %in% c('INS', 'GLOB', 'VASP', 'RRAS', 'TNFRSF12A', 'GPSM3', 'EPHA2', 'PPL', 'EPHA2')) %>%
  mutate(Sample = paste(Sample, ' ChIP')) %>%
  mutate(Line = ifelse(Line == 'sh90', 'shARID1A-2', as.character(Line))) %>%
  mutate(Line = ifelse(Line == 'sh91', 'shARID1A-1', as.character(Line))) %>%
  mutate(Line = factor(Line, levels = c('shNS', 'shARID1A-1', 'shARID1A-2'))) %>%
  ggplot(aes(x = Primer, y = u, fill = Line)) + 
  geom_errorbar(aes(ymin = u-se, ymax =u +se), position = position_dodge(0.9), width = 0.2) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  geom_bar(position = 'dodge', stat = 'identity', color = 'grey30', show_guide = F) + 
  scale_fill_manual(values = c('grey50', 'steelblue', 'steelblue3')) + 
  facet_wrap(~Sample, scales = 'free_y', ncol = 1) + theme_paper() + 
  xlab('') + ylab('Fold Enrichment Relative to INS Promoter') + 
  scale_y_continuous(expand = c(0,.1)) +
  geom_hline(y = 0, color = 'grey10') +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_line(color = 'grey80', linetype = 'dashed'), 
        axis.text = element_text(color = 'grey10'),
        axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        strip.background = element_rect(colour = "grey20", fill = NA)  )
ggsave('~/proj/raab_swisnf_2015/figures/chip.pdf', plt, width = 10, height = 10)

# T.test on samples with multiple


sh90_pvals <- comb %>%
  filter(!Primer %in% c('INS')) %>%
  filter(Sample %in% c('H3K27ac', 'ARID1A', 'ARID2'), 
         Primer %in% c('VASP', 'RRAS', 'TNFRSF12A', 'EPHA2', 'GPSM3', 'NKD1', 'PPL', 'EPHA2'), 
         Line %in% c('shNS', 'sh90')) %>%
  ungroup() %>%
  group_by(Sample, Primer) %>%
  do(tidy(t.test(norm_percent ~ Line, data = .) ) )
View(sh90_pvals) 

sh91_pvals <- comb %>%
  filter(!Primer %in% c('INS')) %>%
  filter(Sample %in% c('H3K27ac', 'ARID1A', 'ARID2'), 
         Primer %in% c('VASP', 'RRAS', 'TNFRSF12A', 'GPSM3', 'EPHA2', 'PPL', 'EPHA2', 'NKD1'), 
         Line %in% c('shNS', 'sh91')) %>%
  ungroup() %>%
  group_by(Sample, Primer) %>%
  do(tidy(t.test(norm_percent ~ Line, data = .) ) ) 
View(sh91_pvals)

comb %>% filter(!Sample %in% c('2', '0.2', '0.02')) %>% group_by(Primer, Sample, exp) %>%
  ungroup() %>% 
  group_by(Line, Primer, Sample) %>%
  summarise(u = mean(norm_percent), se = sd(norm_percent)/sqrt(n())) %>%
  filter(Sample %in% c('ARID1A', 'ARID2', 'H3K27ac'), 
         Primer %in% c('INS', 'GLOB', 'VASP', 'RRAS', 'TNFRSF12A', 'GPSM3', 'EPHA2', 'PPL', 'EPHA2', 'NKD1')) %>%
  ggplot(aes(x = Primer, y = u, fill = Line)) + 
  geom_errorbar(aes(ymin = u-se, ymax =u +se), position = position_dodge(0.9), width = 0.2) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  geom_bar(position = 'dodge', stat = 'identity', color = 'grey30', show_guide = F) + 
  scale_fill_manual(values = c('steelblue', 'steelblue3', 'grey50')) + 
  facet_wrap(~Sample, scales = 'free_y', nrow = 1) + theme_paper() + 
  xlab('') + ylab('Fold Enrichment Relative to INS Promoter') + 
  scale_y_continuous(expand = c(0,.1)) +
  geom_hline(y = 0, color = 'grey10') +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_line(color = 'grey80', linetype = 'dashed'), 
        axis.text = element_text(color = 'grey10'),
        axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        strip.background = element_rect(colour = "grey20", fill = NA)  )
