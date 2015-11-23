# Figure 10 
library(ggplot2) 
library(scales) 
library(data.table) 
library(RColorBrewer) 
library(dplyr) 
library(broom)
library(parallel)
library(tidyr)
library(gridExtra)
source('code/util/theme_paper.R')
source('code/misc/cleanNames.R')
source('code/function_defs/fig9_helpers_func.R')

#setup directories
outdir_figs <- 'figures/'
all_h5 <- list.files('output/encode_coverages/coverages/') # requires all .h5 files

################################################################################
#  Look at the directly regulated genes alone
#  Determine if Arid1b and Arid2 directly regulated genes 
#  Separated by Direction and/or genomic location 
#  Have any interesting of differential enrichments 
#
#  Similar analysis for Arid1a and Arid2 except look at 
#  Competitive vs Cooperative Interactions 
# 
#  Global analysis to determine whats intersting 
#  Try metaplots to visualize a few of the best 
#  Make a statement about the regulation of these sets of genes. 
#
################################################################################
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
arid1a_changed <- read.csv('output/diffexp/tables/arid1a_sig_results.csv') %>% add_rownames()
arid1b_changed <- read.csv('output/diffexp/tables/arid1b_sig_results.csv') %>% add_rownames()
arid2_changed <- read.csv('output/diffexp/tables/arid2_sig_results.csv') %>% add_rownames()

# ----------------------------------------------------------------------------------
# Do the Arid1b analysis first 

# Arid1b peaks associated with repressed genes 
# genes repressed by both 
both_repr <- intersect(arid1b_changed[arid1b_changed$log2FoldChange > 0, ]$rowname, arid2_changed[arid2_changed$log2FoldChange > 0 ,]$rowname) 
#506
# genes activated by both 
both_act <- intersect(arid1b_changed[arid1b_changed$log2FoldChange < 0, ]$rowname, arid2_changed[arid2_changed$log2FoldChange < 0, ]$rowname) 
#170

# Find the peaks associated with direct regulation for each mode of interaction
arid1b_repr_direct <- arid1b_peaks %>% 
  filter(nearestTSS %in% both_repr)  

arid2_repr_direct <- arid2_peaks %>%
  filter(nearestTSS %in% both_repr)  

arid1b_act_direct <- arid1b_peaks %>%
  filter(nearestTSS %in% both_act ) 

arid2_act_direct <- arid2_peaks %>% 
  filter(nearestTSS %in% both_act) 



arid1b_2_namelist <- list(arid1b_act_direct$name, arid1b_repr_direct$name, arid2_act_direct$name, arid2_repr_direct$name) 


# Filter the full data frome frome above. 
arid1b_rdv <- tbl_df(arid1b_repr_direct) %>% 
  mutate(group = 'ARID1B', dir = 'Repressed') 
arid2_rdv <- tbl_df(arid2_repr_direct) %>%
  mutate(group = 'ARID2', dir = 'Repressed') 

arid1b_rda <- tbl_df(arid1b_act_direct) %>% 
  mutate(group = 'ARID1B', dir ='Activated') 
arid2_rda <- tbl_df(arid2_act_direct) %>% 
  mutate(group = 'ARID2', dir = 'Activated') 

# Combine data and get correct type information
all_1b_2 <- rbind(arid1b_rdv, arid2_rdv, arid1b_rda, arid2_rda) 
all_1b_2_total <- all_1b_2 %>% 
  mutate(group2 = ifelse(type %in% c('Intron', 'Distal'), 'Enhancer', 'Genic'))

# Do the whole thing over with Arid1a
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

# Filter out the direct targets using this list
arid1a_coop_dp <- arid1a_peaks %>% 
  filter(nearestTSS %in% coop$rowname) 
arid1a_comp_dp <- arid1a_peaks %>%
  filter(nearestTSS %in% comp$rowname)
arid2_coop_dp <- arid2_peaks %>%
  filter(nearestTSS %in% coop$rowname) 
arid2_comp_dp <- arid2_peaks %>%
  filter(nearestTSS %in% comp$rowname)

arid1a_2_namelist <- list(arid1a_coop_dp$name, arid1a_comp_dp$name, 
                          arid2_coop_dp$name, arid2_comp_dp$name) 

arid1a_coop_d <- tbl_df(arid1a_coop_dp) %>% 
  mutate(group = 'ARID1A', dir = 'Cooperative') 
arid2_coop_d <- tbl_df(arid2_coop_dp) %>%
  mutate(group = 'ARID2', dir = 'Cooperative') 

arid1a_comp_d <- tbl_df(arid1a_comp_dp) %>% 
  mutate(group = 'ARID1A', dir ='Competitive') 
arid2_comp_d <- tbl_df(arid2_comp_dp) %>% 
  mutate(group = 'ARID2', dir = 'Competitive') 

all_1a_2 <- rbind(arid1a_coop_d, arid2_coop_d, arid1a_comp_d, arid2_comp_d) 
all_1a_2_total <- all_1a_2 %>%
  mutate(group2 = ifelse(type %in% c('Distal', 'Intron'), 'Enhancer', 'Genic') )

# Calculate pvals for all groups - compare active vs repressed or coop vs competetive
# -----------------------------------------------------------------------------
# Can I take a list of all coverage files and calculate the pvals accordingly. 
all_tails <- unique(unlist(lapply(all_h5, function(x) gsub(pattern = '.*_p_', '', x ))) )
stat_test_cvc <- lapply(all_tails, function(x, nlist) { 
  name = unlist(strsplit(x, split ='_'))[1]
  row = df_for_mp_1a_2(x, nlist) %>% make_bxplot(retPlot = F) %>% ungroup()
  row$name = name   
  return(row)
}, nlist = arid1a_2_namelist)

stat_test_avr <- lapply(all_tails, function(x, nlist) { 
  name = unlist(strsplit(x, split ='_'))[1]
  row = df_for_mp_1b_2(x, nlist) %>% make_bxplot(retPlot = F) %>% ungroup()
  row$name = name   
  return(row)
}, nlist = arid1b_2_namelist)

enrich_cvc <- lapply(all_tails, function(x, nlist) { 
  name = unlist(strsplit(x, split = '_'))[1]
  row = df_for_mp_1a_2(x, nlist) %>% make_signal_table() %>% ungroup()
  row$name =  name
  return(row)
}, nlist = arid1a_2_namelist )

enrich_avr <- lapply(all_tails, function(x, nlist) { 
  name = unlist(strsplit(x, split = '_'))[1]
  row = df_for_mp_1b_2(x, nlist) %>% make_signal_table() %>% ungroup()
  row$name =  name
  return(row)
}, nlist = arid1b_2_namelist )

cvc_all <- rbind_all(stat_test_cvc) %>% mutate(padj = p.adjust(p.value, 'BH') ) 
avr_all <- rbind_all(stat_test_avr)  %>% mutate(padj = p.adjust(p.value, 'BH') )
cvc_sig <- rbind_all(enrich_cvc)
cvc_sig$name <- fixNames(cvc_sig$name) 
cvc_sig <- cvc_sig %>% filter(!name %in% filterNames)  
avr_sig <- rbind_all(enrich_avr) 
avr_sig$name <- fixNames(avr_sig$name) 
avr_sig <- avr_sig %>% filter(!name %in% filterNames) 

# This prints the list of most sig changed
cvc_all %>% filter(padj < 0.05) %>% arrange(padj) 
avr_all %>% filter(padj < 0.05)  %>% arrange(padj) 

# Create plots for Supplemetnal Figure 10-1 and 10-2
cvc_plot <- cvc_sig %>% make_lineplot()
ggsave('figures/supplemental/S9_2.pdf', cvc_plot, scale =2)
avr_plot <- avr_sig %>% make_lineplot()
ggsave('figures/supplemental/S9_1.pdf', avr_plot, scale =2)

# Present most interesting factors in main Figure 10
# ------------------------------------------------------------------------------
brca1_avr_meta <- df_for_mp_1b_2('Brca1a300IggrabAln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() 
brca1_avr_bx <- df_for_mp_1b_2('Brca1a300IggrabAln_s_coverage.h5', arid1b_2_namelist) %>% make_bxplot(retPlot=T) 
avr_all %>% filter(name == 'Brca1a300IggrabAln')

max_avr_meta <- df_for_mp_1b_2('MaxIggrabAln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() +ylab('')
max_avr_bx <- df_for_mp_1b_2('MaxIggrabAln_s_coverage.h5', arid1b_2_namelist) %>% make_bxplot(retPlot=T) +ylab('')
avr_all %>% filter(name == 'MaxIggrabAln')

tead_avr_meta <- df_for_mp_1b_2('Tead4sc101184V0422111Aln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() +ylab('')#
tead_avr_bx <- df_for_mp_1b_2('Tead4sc101184V0422111Aln_s_coverage.h5', arid1b_2_namelist) %>% make_bxplot(retPlot=T) +ylab('')
avr_all %>% filter(name == 'Tead4sc101184V0422111Aln')

tcf_avr_meta <- df_for_mp_1b_2('Tcf7l2UcdAln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() +ylab('')
tcf_avr_bx <- df_for_mp_1b_2('Tcf7l2UcdAln_s_coverage.h5', arid1b_2_namelist) %>% make_bxplot(retPlot=T) +ylab('')
avr_all %>% filter(name == 'Tcf7l2UcdAln')

pdf('figures/Fig_9A_arid1b_arid2_factors.pdf', height = 12, width = 24)
grid.arrange(arrangeGrob(brca1_avr_meta, 
                         max_avr_meta, 
                         tead_avr_meta, 
                         tcf_avr_meta, nrow =1), 
             arrangeGrob(brca1_avr_bx,
                         max_avr_bx, 
                         tead_avr_bx,
                         tcf_avr_bx, nrow = 1), 
             nrow = 2) 
dev.off()

# Do a similar thing using the Arid1a and Arid2 set (competition vs cooperation)
fosl_cvc_meta <- df_for_mp_1a_2('Fosl2V0416101Aln_s_coverage.h5', arid1a_2_namelist) %>% make_meta()
fosl_cvc_bx <- df_for_mp_1a_2('Fosl2V0416101Aln_s_coverage.h5', arid1a_2_namelist) %>% make_bxplot()
cvc_all %>% filter(name == 'Fosl2V0416101Aln')

corest_cvc_meta <- df_for_mp_1a_2('Corestsc30189IggrabAln_s_coverage.h5', arid1a_2_namelist) %>% make_meta()
corest_cvc_bx <- df_for_mp_1a_2('Corestsc30189IggrabAln_s_coverage.h5', arid1a_2_namelist) %>% make_bxplot()
cvc_all %>% filter(name == 'Corestsc30189IggrabAln')

maz_cvc_meta <- df_for_mp_1a_2('Mazab85725IggrabAln_s_coverage.h5', arid1a_2_namelist) %>% make_meta()
maz_cvc_bx <- df_for_mp_1a_2('Mazab85725IggrabAln_s_coverage.h5', arid1a_2_namelist) %>% make_bxplot()
cvc_all %>% filter(name == 'Mazab85725IggrabAln')

cjun_cvc_meta <- df_for_mp_1a_2('CjunIggrabAln_s_coverage.h5', arid1a_2_namelist) %>% make_bxplot(notch=T)

pdf('figures/Fig_9B_arid1a_arid2_factors.pdf', height = 12, width = 18)
grid.arrange(arrangeGrob(fosl_cvc_meta, corest_cvc_meta, maz_cvc_meta, nrow = 1),  
             arrangeGrob(fosl_cvc_bx, corest_cvc_bx, maz_cvc_bx, nrow = 1), nrow = 2)
dev.off()
