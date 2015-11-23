#compare histone enrichment data at peaks
library(ggplot2) 
library(scales) 
library(data.table) 
library(RColorBrewer) 
library(dplyr) 
library(broom)
source('code/theme_paper.R')
#setup directories
HOME <- '~/proj/swi_snf_final/' #should be project home
outdir_figs <- paste0(HOME, 'output/plots/')
outdir_tabs <- paste0(HOME, 'output/plots')

#get enrichment information
edir='output/encode_coverages/enrichments/'
files = list.files(edir)
groups = c('arid1a', 'arid1b', 'arid2')
df = data.frame(gene=NULL, enrichment=NULL, peakname=NULL, signame=NULL)
for (g in groups){ 
  for (f in files) { 
    peak = unlist(strsplit(f, '_'))[1]
    if (peak %in% g) { 
      dat = fread(paste0(edir,f) )
      sig = unlist(strsplit(f, '_'))[3]
      dat$peakname = rep(peak, nrow(dat))
      dat$signame  = rep(sig,  nrow(dat))
      dat$enrichment= as.numeric(as.character(dat$enrichment))
      df = rbind(df, dat)
    }
  }
}
setnames(df, 'gene', 'name')

#fix the names of the antibodies
#several have other annoying info these have to be dealt with one at a time
df$signame = gsub(pattern = 'Igg.*$|Ucd.*$|Std.*$|Aln.*$', replacement = '', x=df$signame, perl = T)
df$signame = gsub(pattern = 'Arid3anb100279', replacement='Arid3a', df$signame, perl=T)
df$signame = gsub(pattern = 'Bhlhe40c', replacement='Bhlhe40', df$signame, perl=T)
df$signame = gsub(pattern = 'Brca1a300', replacement='Brca1', df$signame, perl=T)
df$signame = gsub(pattern = 'Chd2ab68301', replacement='Chd2', df$signame, perl=T)
df$signame = gsub(pattern = 'Corestsc30189', replacement='Corest', df$signame, perl=T)
df$signame = gsub(pattern = 'Ezh239875', replacement='Ezh2', df$signame, perl=T)
df$signame = gsub(pattern = 'Maffm8194', replacement='Maff', df$signame, perl=T)
df$signame = gsub(pattern = 'Mafkab50322', replacement='Mafk', df$signame, perl=T)
df$signame = gsub(pattern = 'Mafksc477', replacement='Mafk_sc', df$signame, perl=T)
df$signame = gsub(pattern = 'Mazab85725', replacement='Maz', df$signame, perl=T)
df$signame = gsub(pattern = 'P300sc582', replacement='P300', df$signame, perl=T)
df$signame = gsub(pattern = 'Rfx5200401194', replacement='Rfx', df$signame, perl=T)
df$signame = gsub(pattern = 'Smc3.*$', replacement = 'Smc3', df$signame, perl=T)
df$signame = gsub(pattern = 'Zeb1.*$', replacement = 'Zeb1', df$signame, perl=T) 
df$signame = gsub(pattern = 'Pcr1x', replacement = '', df$signame, perl=T) 
df$signame = gsub(pattern = 'Pcr2x', replacement = '', df$signame, perl=T) 
df$signame = gsub(pattern = 'Yy1.*$', replacement = 'Yy1', df$signame, perl=T) 
df$signame = gsub(pattern = 'Atf3.*$', replacement = 'Atf3', df$signame, perl=T) 
df$signame = gsub(pattern = 'Bhlhe40V.*$', replacement = 'Bhlhe40_h', df$signame, perl=T)
df$signame = gsub(pattern = 'Cebpd.*$', replacement = 'Cebpd', df$signame, perl=T) 
df$signame = gsub(pattern = 'Creb1.*$', replacement = 'Creb1', df$signame, perl = T) 
df$signame = gsub(pattern = 'Ctcfsc5916V.*$', replacement = 'Ctcf_h', df$signame, perl = T) 
df$signame = gsub(pattern = 'Elf1.*$', replacement = 'Elf1', df$signame, perl=T) 
df$signame = gsub(pattern = 'Fosl2.*$', replacement = 'Fosl2', df$signame, perl = T) 
df$signame = gsub(pattern = 'Foxa1sc101.*$', replacement ='Foxa1_a', df$signame, perl = T) 
df$signame = gsub(pattern = 'Foxa1sc655.*$', replacement = 'Foxa1_b', df$signame, perl = T) 
df$signame = gsub(pattern = 'Hdac2.*$', replacement = 'Hdac2', df$signame, perl =T) 
df$signame = gsub(pattern = 'Hey1.*$', replacement = 'Hey1', df$signame, perl=T) 
df$signame = gsub(pattern = 'Hnf4a.*$', replacement = 'Hnf4a', df$signame, perl = T) 
df$signame = gsub(pattern = 'Hnf4g.*$', replacement = 'Hnf4g', df$signame, perl = T) 
df$signame = gsub(pattern = 'MaxV.*$', replacement = 'Max_h', df$signame, perl = T) 
df$signame = gsub(pattern = 'Mbd4.*$', replacement = 'Mbd4', df$signame, perl=T) 
df$signame = gsub(pattern = 'Mybl2.*$', replacement = 'Mybl2', df$signame, perl = T) 
df$signame = gsub(pattern = 'Nfic.*$', replacement = 'Nfic', df$signame, perl = T) 
df$signame = gsub(pattern = 'Nr2f2.*$', replacement = 'Nr2f2', df$signame, perl = T) 
df$signame = gsub(pattern = 'NrsfV.*$', replacement = 'Nrsf_h', df$signame, perl = T) 
df$signame = gsub(pattern = 'P300V.*$', replacement = 'P300_h', df$signame, perl = T) 
df$signame = gsub(pattern = 'Rad21V.*$', replacement = 'Rad21_h', df$signame, perl = T) 
df$signame = gsub(pattern = 'Sin3a.*$', replacement = 'Sin3a', df$signame, perl = T) 
df$signame = gsub(pattern = 'Sp2.*$', replacement = 'Sp2', df$signame, perl = T) 
df$signame = gsub(pattern = 'Srf.*$', replacement = 'Srf1', df$signame, perl = T) 
df$signame = gsub(pattern = 'Tead4.*$', replacement = 'Tead4', df$signame, perl =T) 


#drop antibodies that had treatedments or the one duplicate (mafk_sc) which both antibodies look similiar in my initial look
dropname <- c('Pol2Forskln', 'Srebp1Insln', 'Pgc1aForskln', 'Mafk_sc', 'Hsf1Forskln', 'Grp20Forskln', 'ErraForskln', 'Rxlch', 'RxlchV0416101', 'RxlchV0422111', 'Bhlhe40_h', 'Nsrf') 
df <- df[!df$signame %in% dropname, ]


#create data frames for each arid
arid1a_df = df[df$peakname == 'arid1a', ]
arid1a_peaks = read.csv('output/macs_peaks/cleaned/arid1a_annotated.csv') 
arid1a_df = merge(arid1a_df, arid1a_peaks, by='name')

arid1b_df = df[df$peakname == 'arid1b', ]
arid1b_peaks = read.csv('output/macs_peaks/cleaned/arid1b_annotated.csv') 
arid1b_df = merge(arid1b_df, arid1b_peaks, by='name')

arid2_df  = df[df$peakname == 'arid2', ]
arid2_peaks = read.csv('output/macs_peaks/cleaned/arid2_annotated.csv') 
arid2_df = merge(arid2_df, arid2_peaks, by='name') 



#combined data frames
full = rbind(arid1a_df, arid1b_df, arid2_df)
full$type <- factor(full$type, levels=c('Distal', 'Intron', 'Promoter', 'Exon') )


full_summary_byclass <- full %>% 
  group_by(peakname, signame, svm) %>%
  summarise(avg=mean(enrichment) )

full_summary_bypos <- full %>% 
  group_by(peakname, signame, type) %>% 
  summarise(avg=mean(enrichment) ) 


#so I can pull out just the histones 
histones = c('H4k20me1', 'H3k9ac', 'H3k79me2', 'H3k4me2', 'H3k4me3', 'H3k36me3', 
             'H3k27me3', 'H3k27ac', 'H3k09me3', 'H3k04me1', 'H2az')

#separate data by histone vs not
histone_class_summary <- full_summary_byclass[full_summary_byclass$signame %in% histones, ]
txn_class_summary <- full_summary_byclass[!full_summary_byclass$signame %in% histones, ] 

histone_type_summary <- full_summary_bypos[full_summary_bypos$signame %in% histones, ] 
txn_type_summary <- full_summary_bypos[!full_summary_bypos$signame %in% histones, ]   

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}

histone_class_summary$peakname <- sapply(histone_class_summary$peakname, simpleCap)
txn_class_summary$peakname <- sapply(txn_class_summary$peakname, simpleCap) 
histone_class_summary$peakname <- sapply(histone_class_summary$peakname, simpleCap)
txn_class_summary$peakname <- sapply(txn_class_summary$peakname, simpleCap) 

makeheat <- function(df, facet) { 
  col= brewer.pal(9,'Reds')
  theme_clean = function() { 
    theme(text=element_text(family='Helvetica'),
          axis.ticks = element_blank(), 
          axis.text.x=element_text(size=16, angle=90, vjust=0.5), 
          axis.text.y=element_text(size=14), 
          legend.title=element_blank(), 
          panel.background=element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(size=16)) 
  }
  p <- ggplot(df, aes(x=peakname, y=signame, fill=avg) ) + geom_tile(colour='grey30') 
  p <- p + facet_wrap(as.formula(facet), nrow=1 ) 
  p <- p + theme_clean() + xlab('') + ylab('') 
  p <- p + scale_fill_gradient2(low='white', high='red')  
  return(p) 
}

hplot_class <- makeheat(histone_class_summary, facet= '~svm' )  +coord_fixed()
hplot_class 
ggsave(paste0(outdir_figs, 'histones_byclass_enrichment.pdf'), hplot_class) 
ggsave(paste0(outdir_figs, 'histones_byclass_enrichment.tiff'), hplot_class, dpi = 200)

tplot_class <- makeheat(txn_class_summary, facet = '~svm') + coord_fixed()
tplot_class 
ggsave(paste0(outdir_figs, 'txnfactors_byclass_enrichment.pdf'), tplot_class) 
ggsave(paste0(outdir_figs, 'txnfactors_byclass_enrichment.tiff'), tplot_class)

hplot_type <- makeheat(histone_type_summary, facet = '~type')
hplot_type 
ggsave(paste0(outdir_figs, 'histones_bytype_enrichment.pdf') , hplot_type)

tplot_type <- makeheat(txn_type_summary, facet = '~type' ) 
tplot_type
ggsave(paste0(outdir_figs, 'txnfactors_bytype_enrichment.pdf'), tplot_type)

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

# What would really by ideal is sorting out peaks associated with activation vs repression 
# Seeing if there are any differences here. 

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
  
arid2_repr_direct <- arid2_peaks %>%
  filter(nearestTSS %in% both_repr)  
  
arid1b_act_direct <- arid1b_peaks %>%
  filter(nearestTSS %in% both_act ) 
 
arid2_act_direct <- arid2_peaks %>% 
  filter(nearestTSS %in% both_act) 
 
# Now I havea  list I can filter the full data frome frome above. 
arid1b_rdv <- tbl_df(df) %>% 
  filter(name %in% arid1b_repr_direct$name)  %>%
  mutate(group = 'Arid1b', dir = 'Repressed') 
arid2_rdv <- tbl_df(df) %>%
  filter(name %in% arid2_repr_direct$name) %>%
  mutate(group = 'Arid2', dir = 'Repressed') 

arid1b_rda <- tbl_df(df) %>% 
  filter(name %in% arid1b_act_direct$name) %>%
  mutate(group = 'Arid1b', dir ='Activated') 
arid2_rda <- tbl_df(df) %>% 
  filter(name %in% arid2_act_direct$name) %>%
  mutate(group = 'Arid2', dir = 'Activated') 



all_1b_2 <- rbind(arid1b_rdv, arid2_rdv, arid1b_rda, arid2_rda) 
all_1b_2_total <- all_1b_2 %>% 
  left_join(peak_lookup, by = 'name') %>%
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

arid1a_coop_d <- tbl_df(df) %>% 
  filter(name %in% arid1a_coop_dp$name) %>%
  mutate(group = 'Arid1a', dir = 'Cooperative') 
arid2_coop_d <- tbl_df(df) %>%
  filter(name %in% arid2_coop_dp$name) %>%
  mutate(group = 'Arid2', dir = 'Cooperative') 

arid1a_comp_d <- tbl_df(df) %>% 
  filter(name %in% arid1a_comp_dp$name) %>%
  mutate(group = 'Arid1a', dir ='Competitive') 
arid2_comp_d <- tbl_df(df) %>% 
  filter(name %in% arid2_comp_dp$name) %>%
  mutate(group = 'Arid2', dir = 'Competitive') 

all_1a_2 <- rbind(arid1a_coop_d, arid2_coop_d, arid1a_comp_d, arid2_comp_d) 
all_1a_2_total <- all_1a_2 %>%
  left_join(peak_lookup, by = 'name') %>%
  mutate(group2 = ifelse(type %in% c('Distal', 'Intron'), 'Enhancer', 'Genic') )


# Metaplots of most interesting groups 
# Metaplots for comparisons between Arid1b and Arid2 directly regulated genes 
# And between ARid1a and ARid2 competitive and cooperative directly regulated genes

# Draws heavily on the original metaplots scripts. 

summarise_data_group <- function(df) { 
  df %>% 
    group_by(signame, peakname, type) %>% 
    summarise(avg = mean(enrichment), med = median(enrichment) ) %>%
    ungroup() %>% 
    arrange(dplyr::desc(avg) ) 
}
# I can use the metaplot hdf5 data directly for testing differences. 

# For these I'll need to subset the peaks and get the signal from the hdf5 files
readh5vals <- function(fn) {
  tmp <- suppressWarnings(h5read(fn, 'coverages', bit64conversion='bit64'))
  vals <- data.frame(t(tmp$block0_values) )
  vals$rowname <- tmp$axis1
  return(vals) 
}

filter_arid_coverage <- function(arid_coverage_file, peak_names) { 
  # will return just the coverage for the correct list of peaks
  vals <- readh5vals(arid_coverage_file) 
  return(vals[vals$rowname %in% peak_names,]) 
}


# Following two functions are needed to create a data.fame 
# with suitable labels 
# In all cases I'm just reading in HDF5 data adding information for which 
# peaks are the differnent classes and smushing them back together 
# I'll end up with two by two comparisons (arid peak type, and direction of regulation) 

# Function to compare active and repressed arid1b and arid2
df_for_mp_1b_2 <- function(end, pname_list) { 
  # relies on a df of all peak info 
  # pass in the  end of the filenam we want and the lists of peak names to get
  # and a length(4) list of 
  arid1b_act <- pname_list[[1]]
  arid1b_rep <- pname_list[[2]]
  arid2_act  <- pname_list[[3]]
  arid2_rep  <- pname_list[[4]]
  arid1b_fn <- paste('output/encode_coverages/coverages/arid1b_p_', end, sep ='') 
  arid2_fn  <- paste('output/encode_coverages/coverages/arid2_p_', end, sep = '') 
  act_1b_vals <- filter_arid_coverage(arid1b_fn, arid1b_act) %>% 
    mutate(group = 'Active', peak = 'Arid1b Peaks')
  rep_1b_vals <- filter_arid_coverage(arid1b_fn, arid1b_rep) %>%
    mutate(group = 'Repressed', peak = 'Arid1b Peaks') 
  act_2_vals  <- filter_arid_coverage(arid2_fn, arid2_act)  %>%
    mutate(group = 'Active', peak = 'Arid2 Peaks') 
  rep_2_vals <- filter_arid_coverage(arid2_fn, arid2_rep) %>%
    mutate(group = 'Repressed', peak = 'Arid2 Peaks') 
  all <- rbind(act_1b_vals, rep_1b_vals, act_2_vals, rep_2_vals) 
  all <- all %>% left_join(peak_lookup, by = c('rowname' = 'name')  )
  all %<>% gather('Pos', 'Value', -group, -rowname, -peak, -type) %>% tbl_df() %>%
    mutate(Pos2 = as.numeric(gsub(pattern = 'X', replacement = '', x = Pos)  ) )
  return(all)
} 


# Function to compare cooperative and competitive 1a and 2
df_for_mp_1a_2 <- function(end, pname_list) { 
  # relies on a df of all peak info 
  # pass in the  end of the filenam we want and the lists of peak names to get
  # and a length(4) list of 
  arid1a_coop <- pname_list[[1]]
  arid1a_comp <- pname_list[[2]]
  arid2_coop  <- pname_list[[3]]
  arid2_comp  <- pname_list[[4]]
  arid1a_fn <- paste('output/encode_coverages/coverages/arid1a_p_', end, sep ='') 
  arid2_fn  <- paste('output/encode_coverages/coverages/arid2_p_', end, sep = '') 
  coop_1a_vals <- filter_arid_coverage(arid1a_fn, arid1a_coop) %>% 
    mutate(group = 'Cooperative', peak = 'Arid1a Peaks')
  comp_1a_vals <- filter_arid_coverage(arid1a_fn, arid1a_comp) %>%
    mutate(group = 'Competitive', peak = 'Arid1a Peaks') 
  coop_2_vals  <- filter_arid_coverage(arid2_fn, arid2_coop)  %>%
    mutate(group = 'Cooperative', peak = 'Arid2 Peaks') 
  comp_2_vals <- filter_arid_coverage(arid2_fn, arid2_comp) %>%
    mutate(group = 'Competitive', peak = 'Arid2 Peaks') 
  all <- rbind(coop_1a_vals, comp_1a_vals, coop_2_vals, comp_2_vals)
  all <- all %>% left_join(peak_lookup, by = c('rowname' = 'name') )
  all %<>% gather('Pos', 'Value', -group, -rowname, -peak, -type) %>% tbl_df() %>%
    mutate(Pos2 = as.numeric(gsub(pattern = 'X', replacement = '', x = Pos)  ) )
  return(all)
} 
# The above two functions work. Can be incorporated into other workflow.

#function to handle the plotting
make_meta <- function(dfr) { 
  # takes the output of the above functions and plots by group (color) and
  # facets by peak type
  # need to clean up the numbering 
  plt <- dfr %>%  
    group_by(peak, group,  Pos2) %>%
    summarise(avg  = mean(Value), ci = 1.96*(sd(Value))/sqrt(n() ) ) %>%
    mutate(ymin = avg - ci, ymax= avg + ci) %>%
    mutate(Position = rep(seq(-2500, 2490, by = 10)) ) %>%
    ggplot(aes(x=Position, y = avg, color = group, fill = group)) + 
    geom_line(size = 1.5) +
    facet_wrap(~peak) + 
    theme_paper() + 
    scale_fill_manual(values = c('grey20', 'red2') )+ 
    scale_color_manual(values  = c('grey20', 'red2'))  + 
    scale_x_continuous(labels = c('-2', '-1', '0', '1', '2'), breaks = c(-2000, -1000, 0, 1000, 2000) ) + 
    theme(panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = NA, color = 'grey30'), 
          legend.position = 'none') + 
    xlab('Distance from Peak Center (Kb)') + ylab('Relative Enrichment') 
  return(plt)
}

make_bxplot <- function(dfr, retPlot) { 
 # just like above but output is a boxplot 
  pltdf <- dfr %>%  
    filter(Pos2 > 150, Pos2 < 350 ) %>%
    group_by(group, peak, rowname)  %>%
    summarise(signal = mean(Value) )
  summary_pvals <- pltdf %>% 
    ungroup() %>%
    group_by(peak) %>% 
    do(tidy(wilcox.test(signal ~ group, data = .)) ) %>%
    print()
  if (retPlot == FALSE){ 
    return(summary_pvals)
    }
  else{  
    plt <- pltdf %>% 
      ggplot(aes(x=group, y = signal, color = group) ) + 
      geom_boxplot(notch =T) +
      facet_wrap(~peak) + 
      theme_paper() + 
      scale_color_manual(values = c('grey20', 'red2') ) + 
      theme(panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = NA, color = 'grey30'), 
          legend.position = 'none', 
          axis.text = element_text(color = 'grey30'), 
          axis.ticks = element_blank() ) + 
      xlab('') + ylab('Median Relative Enrichment') 
    
    return(plt)
  }
}

arid1b_2_namelist <- list(arid1b_act_direct$name, arid1b_repr_direct$name, arid2_act_direct$name, arid2_repr_direct$name) 
# numbers 
sapply(arid1b_2_namelist, length) 
# 88, #287, #119, #489

arid1a_2_namelist <- list(arid1a_coop_dp$name, arid1a_comp_dp$name, 
                          arid2_coop_dp$name, arid2_comp_dp$name) 
sapply(arid1a_2_namelist, length) 
# 63, 37, 69 , 23
# This works well
# Now all I need to do is put in the stubs and make and save comparisons


# Can I take a list of all coverage files and calculate the pvals accordingly. 
all_h5 <- list.files('output/encode_coverages/coverages/')
all_tails <- unique(unlist(lapply(all_h5, function(x) gsub(pattern = '.*_p_', '', x ))) )

stat_test_cvc <- lapply(all_tails, function(x, nlist) { 
    name = unlist(strsplit(x, split ='_'))[1]
    row = df_for_mp_1a_2(x, nlist) %>% make_bxplot(retPlot = F) %>% ungroup()
    row$name = ame   
    return(row)
  }, nlist = arid1a_2_namelist)

stat_test_avr <- lapply(all_tails, function(x, nlist) { 
  name = unlist(strsplit(x, split ='_'))[1]
  row = df_for_mp_1b_2(x, nlist) %>% make_bxplot(retPlot = F) %>% ungroup()
  row$name = name   
  return(row)
}, nlist = arid1b_2_namelist)

cvc_all <- rbind_all(stat_test_cvc) %>% mutate(padj = p.adjust(p.value, 'BH') )
avr_all <- rbind_all(stat_test_avr)  %>% mutate(padj = p.adjust(p.value, 'BH') )


cvc_all %>% filter(padj < 0.05)
avr_all %>% filter(padj < 0.05)

df_for_mp_1b_2('Srebp1InslnStdAln_s_coverage.h5', arid1b_2_namelist) %>%  make_meta() #
df_for_mp_1b_2('Srebp1InslnStdAln_s_coverage.h5', arid1b_2_namelist) %>%  make_bxplot(retPlot =T) #

df_for_mp_1b_2('Usf2IggrabAln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() 
df_for_mp_1b_2('Usf2IggrabAln_s_coverage.h5', arid1b_2_namelist) %>% make_bxplot(retPlot=T) 



df_for_mp_1b_2('Tead4sc101184V0422111Aln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() #
df_for_mp_1b_2('Usf1Pcr1xAln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() 
df_for_mp_1b_2('Sp1Pcr1xAln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() 
df_for_mp_1b_2('RxraPcr1xAln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() 
df_for_mp_1b_2('H3k27acStdAln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() 
df_for_mp_1b_2('H3k27me3StdAln_s_coverage.h5', arid1b_2_namelist) %>% make_meta()
df_for_mp_1b_2('Foxa1sc6553V0416101Aln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() 
df_for_mp_1b_2('H3k04me1StdAln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() 
df_for_mp_1b_2('H3k4me3StdAln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() 
df_for_mp_1b_2('Mybl2sc81192V0422111Aln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() #
df_for_mp_1b_2('Nficsc81335V0422111Aln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() #
df_for_mp_1b_2('Mbd4sc271530V0422111Aln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() 
df_for_mp_1b_2('Hnf4asc8987V0416101Aln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() #
df_for_mp_1b_2('Hnf4gsc6558V0416101Aln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() #
df_for_mp_1b_2('MaxIggrabAln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() 
df_for_mp_1b_2('Yy1sc281V0416101Aln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() 
df_for_mp_1b_2('Arid3anb100279IggrabAln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() 
df_for_mp_1b_2('Sin3ak20Pcr1xAln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() 
df_for_mp_1b_2('Bhlhe40V0416101Aln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() 
df_for_mp_1b_2('Bhlhe40cIggrabAln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() 
df_for_mp_1b_2('Brca1a300IggrabAln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() 
df_for_mp_1b_2('Corestsc30189IggrabAln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() 
df_for_mp_1b_2('CebpbIggrabAln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() 
df_for_mp_1b_2('Corestsc30189IggrabAln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() 
df_for_mp_1b_2('Cebpdsc636V0416101Aln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() 
df_for_mp_1b_2('CebpzIggrabAln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() 
df_for_mp_1b_2('Ezh239875Aln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() 
df_for_mp_1b_2('H2azStdAln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() 
df_for_mp_1b_2('JundIggrabAln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() 
df_for_mp_1b_2('Tr4UcdAln_s_coverage.h5', arid1b_2_namelist) %>% make_meta() 

df_for_mp_1b_2('Tr4UcdAln_s_coverage.h5', arid1b_2_namelist) %>% make_bxplot() 
df_for_mp_1b_2('JundIggrabAln_s_coverage.h5', arid1b_2_namelist) %>% make_bxplot() 



df_for_mp_1a_2('CtcfStdAln_s_coverage.h5', arid1a_2_namelist) %>% make_meta() 
df_for_mp_1a_2('MaxIggrabAln_s_coverage.h5', arid1a_2_namelist) %>% make_meta()
df_for_mp_1a_2('Yy1sc281V0416101Aln_s_coverage.h5', arid1a_2_namelist) %>% make_meta()
df_for_mp_1a_2('Arid3anb100279IggrabAln_s_coverage.h5', arid1a_2_namelist) %>% make_meta()
df_for_mp_1a_2('Foxa1sc6553V0416101Aln_s_coverage.h5', arid1a_2_namelist) %>% make_meta()
df_for_mp_1a_2('RxraPcr1xAln_s_coverage.h5', arid1a_2_namelist) %>% make_meta()
df_for_mp_1a_2('H3k27acStdAln_s_coverage.h5', arid1a_2_namelist) %>% make_meta()
df_for_mp_1a_2('H3k04me1StdAln_s_coverage.h5', arid1a_2_namelist) %>% make_meta()
df_for_mp_1a_2('H3k4me3StdAln_s_coverage.h5', arid1a_2_namelist) %>% make_meta()
df_for_mp_1a_2('Sp1Pcr1xAln_s_coverage.h5', arid1a_2_namelist) %>% make_meta()
df_for_mp_1a_2('Mazab85725IggrabAln_s_coverage.h5', arid1a_2_namelist) %>% make_meta()
df_for_mp_1a_2('CjunIggrabAln_s_coverage.h5', arid1a_2_namelist) %>% make_meta()
df_for_mp_1a_2('Rfx5200401194IggrabAln_s_coverage.h5', arid1a_2_namelist) %>% make_meta()
df_for_mp_1a_2('H3k27me3StdAln_s_coverage.h5', arid1a_2_namelist) %>% make_meta()
df_for_mp_1a_2('Ezh239875Aln_s_coverage.h5', arid1a_2_namelist) %>% make_meta()
df_for_mp_1a_2('H2azStdAln_s_coverage.h5', arid1a_2_namelist) %>% make_meta()
df_for_mp_1a_2('Fosl2V0416101Aln_s_coverage.h5', arid1a_2_namelist) %>% make_meta()
df_for_mp_1a_2('JundIggrabAln_s_coverage.h5', arid1a_2_namelist) %>% make_meta()
df_for_mp_1a_2('JundPcr1xAln_s_coverage.h5', arid1a_2_namelist) %>% make_meta()
df_for_mp_1a_2('Fosl2V0416101Aln_s_coverage.h5', arid1a_2_namelist) %>% make_meta()
df_for_mp_1a_2('Irf3IggrabAln_s_coverage.h5', arid1a_2_namelist) %>% make_meta()




df_for_mp_1a_2('H3k27acStdAln_s_coverage.h5', arid1a_2_namelist) %>% make_bxplot()
df_for_mp_1a_2('H3k04me1StdAln_s_coverage.h5', arid1a_2_namelist) %>% make_bxplot()
df_for_mp_1a_2('H3k4me3StdAln_s_coverage.h5', arid1a_2_namelist) %>% make_bxplot()
df_for_mp_1a_2('Foxa1sc6553V0416101Aln_s_coverage.h5', arid1a_2_namelist) %>%make_bxplot() + theme(legend.position = 'top') 
df_for_mp_1a_2('Irf3IggrabAln_s_coverage.h5', arid1a_2_namelist) %>% make_bxplot()
df_for_mp_1a_2('Fosl2V0416101Aln_s_coverage.h5', arid1a_2_namelist) %>% make_bxplot()


arid1b_summ_act <- arid1b_peaks %>% filter(name %in% arid1b_2_namelist[[1]]) %>% 
  group_by(type) %>% summarise(count = n() ) %>% 
  mutate(percent = count/sum(count), group = 'Arid1b_Activated' )

arid1b_summ_rep <- arid1b_peaks %>% filter(name %in% arid1b_2_namelist[[2]]) %>%
  group_by(type) %>% summarise(count = n() ) %>%
  mutate(percent = count/sum(count), group = 'Arid1b_Repressed' )

arid2_summ_act <- arid2_peaks %>% filter(name %in% arid1b_2_namelist[[3]]) %>% 
  group_by(type) %>% summarise(count = n() ) %>% 
  mutate(percent = count/sum(count), group = 'Arid2_Activated') 

arid2_summ_rep <- arid2_peaks %>% filter(name %in% arid1b_2_namelist[[4]]) %>% 
  group_by(type) %>% summarise(count = n() ) %>% 
  mutate(percent = count/sum(count), group = 'Arid2_Repressed') 
all_1b_two_peaks <- rbind(arid1b_summ_act, arid1b_summ_rep, arid2_summ_act, arid2_summ_rep) %>%
  separate(group, into = c('group', 'dir'), sep = '_') 
all_1b_two_peaks %>% 
  ggplot(aes(x= dir, y = percent, fill = type) ) +
  geom_bar(position = 'stack', stat = 'identity')  + 
  facet_wrap(~group, nrow = 1) 

all_1b_two_peaks %>% filter(group == 'Arid1b') %>% select(count, dir, type) %>% spread(dir, count) %>% select(Activated, Repressed) %>% chisq.test()
# 1b is significant 0.001
all_1b_two_peaks %>% filter(group == 'Arid2') %>% select(count, dir, type) %>% spread(dir, count) %>% select(Activated, Repressed) %>% chisq.test()
# Arid2 is signifcant at 0.015

arid1a_summ_coop <- arid1a_peaks %>% filter(name %in% arid1a_2_namelist[[1]]) %>%
  group_by(type) %>% summarise(count = n() ) %>%
  mutate(percent = count/sum(count), group = 'Arid1a_Cooperative') 
arid1a_summ_comp <- arid1a_peaks %>% filter(name %in% arid1a_2_namelist[[2]]) %>%
  group_by(type) %>% summarise(count = n() ) %>%
  mutate(percent  = count/sum(count), group = 'Arid1a_Competitive') 
arid2_summ_coop <- arid2_peaks %>% filter(name %in% arid1a_2_namelist[[3]]) %>%
  group_by(type) %>% summarise(count = n() ) %>%
  mutate(percent = count/sum(count), group = 'Arid2_Cooperative') 
arid2_summ_comp <- arid2_peaks %>% filter(name %in% arid1a_2_namelist[[4]]) %>%
  group_by(type) %>% summarise(count = n() ) %>%
  mutate(percent = count/sum(count), group = 'Arid2_Competitive') 

all_1a_2_peaks <- rbind(arid1a_summ_coop, arid1a_summ_comp, arid2_summ_coop, arid2_summ_comp) %>%
  separate(group, into = c('group', 'dir'), sep = '_') %>%
  mutate(type = factor(type, levels = c('Promoter', 'Exon', 'Distal', 'Intron') ) )

all_1a_2_peaks %>%
  ggplot(aes(x= dir, y = percent, fill = type) ) +
  geom_bar(position = 'stack', stat = 'identity')  + 
  facet_wrap(~group, nrow = 1)
all_1a_2_peaks %>% filter(group == 'Arid1a') %>% select(count, dir, type) %>% spread(dir, count) %>% select(Competitive, Cooperative) %>% chisq.test()
# no diff for arid1a
all_1a_2_peaks %>% filter(group == 'Arid2') %>% select(count, dir, type) %>% spread(dir, count) %>% select(Competitive, Cooperative) %>% chisq.test()
# no diff for arid2
