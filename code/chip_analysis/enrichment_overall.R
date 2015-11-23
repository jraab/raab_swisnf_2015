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