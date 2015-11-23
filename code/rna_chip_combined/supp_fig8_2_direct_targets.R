# What is the overlap of peaks associated with directly regulated genes. 
# Does it vary based on geonmic feature or type of regulation 

library(dplyr) 
library(ggplot2)
library(tidyr) 
library(broom) 
source('code/theme_paper.R')
source('~/helper_scripts/code/chippeak_annotation.R')

# What would really by ideal is sorting out peaks associated with activation vs repression 
# Seeing if there are any differences here. 
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}
# Raw data I'll need

# Peaks
arid1a_peaks <- read.csv('output/macs_peaks/cleaned/arid1a_annotated.csv')
arid1b_peaks <- read.csv('output/macs_peaks/cleaned/arid1b_annotated.csv')
arid2_peaks  <- read.csv('output/macs_peaks/cleaned/arid2_annotated.csv')
multi_peaks  <- read.csv('output/macs_peaks/cleaned/multi_bound_peaks.csv', header =F)
colnames(multi_peaks) <- c('chr', 'start', 'end', 'name', 'nearestTSS', 'distance', 'state')
all_peaks <- rbind(arid1a_peaks, arid1b_peaks, arid2_peaks) 
peak_lookup <- all_peaks %>% select(name, type) 

# Expression  
arid1a_all <- read.csv('output/diffexp/tables/arid1a_full_results.csv') %>% add_rownames()
arid1a_changed <- read.csv('output/diffexp/tables/arid1a_sig_results.csv') %>% add_rownames()
arid1b_all <- read.csv('output/diffexp/tables/arid1b_full_results.csv') %>% add_rownames()
arid1b_changed <- read.csv('output/diffexp/tables/arid1b_sig_results.csv') %>% add_rownames()
arid2_all  <- read.csv('output/diffexp/tables/arid2_full_results.csv') %>% add_rownames()
arid2_changed <- read.csv('output/diffexp/tables/arid2_sig_results.csv') %>% add_rownames()

# Do the Arid1b analysis first 
# Arid1b peaks associated with repressed genes 
# genes repressed by both 
both_repr <- intersect(arid1b_changed[arid1b_changed$log2FoldChange > 0, ]$rowname, arid2_changed[arid2_changed$log2FoldChange > 0 ,]$rowname) 
#506
# genes activated by both 
both_act <- intersect(arid1b_changed[arid1b_changed$log2FoldChange < 0, ]$rowname, arid2_changed[arid2_changed$log2FoldChange < 0, ]$rowname) 
#170

# Now filter out the direct targets using this list

arid1b_repr_direct <- arid1b_peaks %>% 
  filter(nearestTSS %in% both_repr)  
 #287 peaks (158 genes)
158/505

arid2_repr_direct <- arid2_peaks %>%
  filter(nearestTSS %in% both_repr)  
 #484 peaks (227 genes)
227/505

arid1b_act_direct <- arid1b_peaks %>%
  filter(nearestTSS %in% both_act ) 
 #88 (58 genes)
arid2_act_direct <- arid2_peaks %>% 
  filter(nearestTSS %in% both_act) 
 #119 (72)
df_to_gr <- function(df, padding) {
  gr <- GRanges(df[,1], IRanges(df[,2]-(padding/2), df[,3] + (padding/2) ) ) 
  elementMetadata(gr) <- df[,4:ncol(df)]
  return(gr) 
} 

length(intersect(df_to_gr(arid1b_repr_direct, 500), df_to_gr(arid2_repr_direct, 500) )) 
# 116
length(intersect(df_to_gr(arid2_repr_direct, 500), df_to_gr(arid1b_repr_direct, 500) ) )
#116 (symettrical)
length(intersect(df_to_gr(arid1b_act_direct, 500), df_to_gr(arid2_act_direct, 500) ) )
#40
length(intersect(df_to_gr(arid2_act_direct, 500), df_to_gr(arid1b_act_direct, 500) ) )
#40 

compare_ov <- function(df1, df2, padding = 500) { 
  gr1 <- df_to_gr(df1, padding) 
  gr2 <- df_to_gr(df2, padding) 
  overlap1 <- gr1[overlapsAny(gr1, gr2),]
  anti_1 <- gr1[!overlapsAny(gr1, gr2), ]
  overlap2 <- gr2[overlapsAny(gr2, gr1), ]
  anti_2 <- gr2[!overlapsAny(gr2,gr1) ]
  ov1_df <- as.data.frame(overlap1) %>%mutate(rel = 'Overlap') 
  anti_1_df <- as.data.frame(anti_1) %>% mutate(rel = 'Not Overlap') 
  ov2_df <- as.data.frame(overlap2) %>% mutate(rel = 'Overlap') 
  anti_2_df <- as.data.frame(anti_2) %>% mutate(rel = 'Not Overlap') 
  return(rbind(ov1_df, ov2_df, anti_1_df, anti_2_df) )
  }

x <- compare_ov(arid1b_repr_direct, arid2_repr_direct) %>% tbl_df() %>% mutate(dir = 'Repressed')

y <- compare_ov(arid1b_act_direct, arid2_act_direct) %>% tbl_df() %>% mutate(dir = 'Activated')

xx <- rbind(x,y)  %>% separate(name, into = c('Arid'), sep = '_', remove = F, extra = 'drop') %>%
  tidyr::expand(Arid, type, dir, rel) %>% 
  left_join(rbind(x,y) ) %>% 
  group_by(Arid, type, dir, rel) %>% 
  summarise(count = n()-1 ) %>% 
  mutate(percent = count/sum(count)) %>% 
  ungroup() %>% 
  mutate(Arid = sapply(Arid, simpleCap) ) 

# Figure 8-2 A
xx %>% 
ggplot(aes(x = type, y = percent, fill = rel)) + 
  geom_bar(position = 'stack', stat = 'identity') + 
  geom_bar(position = 'stack', stat = 'identity', color = 'grey30', show_guide = FALSE) +
  facet_grid(Arid ~ dir) + 
  theme_paper() + 
  theme(axis.ticks= element_blank(), 
        axis.text.x = element_text(angle = 30, hjust = 1, vjust =1, color = 'grey20' )) +
  scale_fill_manual(values = c('grey40', 'steelblue')) + 
  xlab('') + ylab('Fraction of Peaks') 

total <- xx %>% 
  group_by(Arid, dir) %>% 
  summarise(count = sum(count) ) %>% 
  ungroup() %>%
  mutate(Arid = 'Total')

print(total)
xx %>% filter(rel == 'Overlap') %>% 
  group_by(Arid, dir) %>% 
  summarise(count = sum(count) ) %>% print() 

overlap <- xx %>% 
  group_by(Arid, dir, rel) %>% 
  summarise(count = sum(count) )
 
# 200 total peaks near activatd genes (50 overlap between 1b and 2 - 25%) 
# 764 total peaks near repressed genes(156 overlap between 1b and 2 - 20.4%) 
overlap %>% filter(Arid == 'Arid1b') %>% 
  select(-Arid) %>% 
  spread(dir, count) %>% print() %>%
  select(-Arid, -rel) %>%
  as.matrix() %>% 
  fisher.test(Activated, Repressed ) %>% 
  tidy()

overlap %>% filter(Arid == 'Arid1b') %>% 
  ggplot(aes(x =dir, y = count, fill = rel)) + 
  geom_bar(position = position_dodge(0.55), stat = 'identity', width = 0.5) + 
  geom_bar(position = position_dodge(0.55), stat = 'identity', color = 'grey30', show_guide =F, width = 0.5) + 
  theme_paper() + 
  scale_fill_manual(values = c('grey40', 'steelblue') ) + 
  scale_y_continuous(expand = c(0,0) )  +
  theme(axis.ticks = element_blank(), 
        panel.grid = element_blank(), 
        axis.line = element_line(), 
        axis.text = element_text(color = 'grey20')) + 
  xlab('') + ylab('Number of Peaks') 

# Signficant difference between proportion active and repressed
# Arid1a and Arid2 

# Plot just half the data becasue it is symmetrical.



arid1a_arid2_direct <- intersect(arid1a_changed$rowname, arid2_changed$rowname) 
arid1a_expr_direct <- arid1a_changed %>% 
  filter(rowname %in% arid1a_arid2_direct) 
arid2_expr_direct <- arid2_changed %>% 
  filter(rowname %in% arid1a_arid2_direct) 
combined_1a_2 <- arid1a_expr_direct %>% 
  inner_join(arid2_expr_direct, by = 'rowname')

comp <- combined_1a_2 %>% 
  filter(sign(log2FoldChange.x) != sign(log2FoldChange.y) )
# 36
coop <- combined_1a_2 %>%
  filter(sign(log2FoldChange.x) == sign(log2FoldChange.y) )
# 59
# not a ton of genes here
# Now I havea  list I can filter the full data frome frome above. 

# First get peaks that are on those lists 
arid1a_coop_dp <- arid1a_peaks %>% 
  filter(nearestTSS %in% coop$rowname) 
arid1a_comp_dp <- arid1a_peaks %>%
  filter(nearestTSS %in% comp$rowname)
arid2_coop_dp <- arid2_peaks %>%
  filter(nearestTSS %in% coop$rowname) 
arid2_comp_dp <- arid2_peaks %>%
  filter(nearestTSS %in% comp$rowname)

x <- compare_ov(arid1a_coop_dp, arid2_coop_dp) %>% tbl_df() %>% mutate(dir = 'Cooperative')
y <- compare_ov(arid1a_comp_dp, arid2_comp_dp) %>% tbl_df() %>% mutate(dir = 'Competetive')

xx <- rbind(x,y)  %>% separate(name, into = c('Arid'), sep = '_', remove = F, extra = 'drop') %>%
  tidyr::expand(Arid, type, dir, rel) %>% 
  left_join(rbind(x,y) ) %>% 
  group_by(Arid, type, dir, rel) %>% 
  summarise(count = n()-1 ) %>% 
  mutate(percent = count/sum(count)) %>% 
  ungroup() %>% 
  mutate(Arid = sapply(Arid, simpleCap) ) 
# 
# Figure 8-2 B
xx %>% 
  ggplot(aes(x = type, y = percent, fill = rel)) + 
  geom_bar(position = 'stack', stat = 'identity') + 
  geom_bar(position = 'stack', stat = 'identity', color = 'grey30', show_guide = FALSE) +
  facet_grid(Arid ~ dir) + 
  theme_paper() + 
  theme(axis.ticks= element_blank(), 
        axis.text.x = element_text(angle = 30, hjust = 1, vjust =1, color = 'grey20' )) +
  scale_fill_manual(values = c('grey40', 'steelblue')) + 
  xlab('') + ylab('Fraction of Peaks') 

total <- xx %>% 
  group_by(Arid, dir) %>% 
  summarise(count = sum(count) ) %>% 
  ungroup() %>%
  mutate(Arid = 'Total')

print(total)
xx %>% filter(rel == 'Overlap') %>% 
  group_by(Arid, dir) %>% 
  summarise(count = sum(count) ) %>% print() 

overlap <- xx %>% 
  group_by(Arid, dir, rel) %>% 
  summarise(count = sum(count) )


overlap %>% filter(Arid == 'Arid1a') %>% 
  select(-Arid) %>% 
  spread(dir, count) %>% print() %>% 
  select(-Arid, -rel) %>% 
  as.matrix() %>% 
  fisher.test(Activated, Repressed ) %>% 
  tidy()

overlap %>% filter(Arid == 'Arid1a') %>% 
  ggplot(aes(x =dir, y = count, fill = rel)) + 
  geom_bar(position = position_dodge(0.55), stat = 'identity', width = 0.5) + 
  geom_bar(position = position_dodge(0.55), stat = 'identity', color = 'grey30', show_guide =F, width = 0.5) + 
  theme_paper() + 
  scale_fill_manual(values = c('grey40', 'steelblue') ) + 
  scale_y_continuous(expand = c(0,0) )  +
  theme(axis.ticks = element_blank(), 
        panel.grid = element_blank(), 
        axis.line = element_line(), 
        axis.text = element_text(color = 'grey20') )  + 
  xlab('') + ylab('Number of Peaks') 
