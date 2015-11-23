# Calculate pair-wise correlation between chip-seq enrichments at arid peaks
# -------------------------------------------------------------------------------
# Script will loop through all enrichment files and calculate the spearman cor 
# between pairs of peaks
# output is saved in in PROJ/output/encode_coverages/ 
# This script takes a few minutes to run so I am separating out from 
#     plot_heatmaps_of_enriched
# --------------------------------------------------------------------------------
library(data.table)
output_dir <- 'output/encode_coverages/' 
enrichment_dir <- 'output/encode_coverages/enrichments/'
source('code/function_defs/cleanNames.R')

# Utility functions  
# -----------------------------------------------------------------------------
getGroupFiles <- function(group, filelist){ 
  filelist_groups <- sapply(filelist, function(x) unlist(strsplit(x, split = '_'))[1] )
  keep <- filelist_groups == group
  print(table(keep) )
  return(filelist[keep]) 
}

getSigname <- function(filename){ 
  return(unlist(strsplit(filename, split = '_'))[3])
}

cleanDF <- function(df) { 
  droppatt <- toupper('Pol2_aForskln|Pol2Forskln*|Srebp1Insln*|Pgc1aForskln*|MAFK_sc*|Hsf1Forskln*|Grp20Forskln*|ErraForskln*|RxlchPcr1x*|RxlchPcr2x*|RxlchV0416101*|RxlchV0422111*')
  df <- df[!grepl(droppatt, x = row.names(df)), !grepl(droppatt, x = colnames(df) )]
  return(df) 
}


# Pairwise r within groups. 
# -----------------------------------------------------------------------------
# Need a list of arid1a, arid1b, arid2 alone peak names 
arid1a_single <- read.csv('output/macs_peaks/cleaned/arid1a_annotated.csv') %>% 
  filter(svm == 'Single') %>% select(name) 
arid1b_single <- read.csv('output/macs_peaks/cleaned/arid1b_annotated.csv') %>% 
  filter(svm == 'Single') %>% select(name) 
arid2_single  <- read.csv('output/macs_peaks/cleaned/arid2_annotated.csv') %>% 
  filter(svm == 'Single') %>% select(name) 
multi_all <- read.csv('output/macs_peaks/cleaned/multi_bound_peaks.csv', header =F) %>% 
  select(V4) 
groups <- c('arid1a', 'arid1b', 'arid2', 'multi')
keeps <- c(arid1a_single, arid1b_single, arid2_single, multi_all)
dropname <- c('Pol2Forskln', 'Srebp1Insln', 'Pgc1aForskln', 'MAFK_sc', 'Hsf1Forskln', 'Grp20Forskln', 'ErraForskln', 'Rxlch', 'RxlchV0416101', 'RxlchV0422111')
files = list.files(enrichment_dir)
output_final <- list(NA, NA, NA, NA)

for (k in 1:length(groups)){  
  g <- groups[k]
  files_interest <- getGroupFiles(g, files) # prints how many files were kept, how many rejected
  # filter out groups int he dropnames 
  files_interest <- files_interest[!files_interest %in% dropname] 
  signal_names <- sapply(files_interest, getSigname) 
  output <- matrix(rep(0), nrow = length(signal_names), ncol = length(signal_names) )
  group_keep <- keeps[[k]]
  # Ideally I would do something more cleaver than reading this data in a million times like below
  #enrichment_list <- lapply(files_interest, function(x) read_in_data(paste0(enrichment_dir, x)) )
  
  # now pairwise through the enrihcment list
  for (f1 in 1:length(files_interest)) {
    for (f2 in 1:length(files_interest)){ 
      dat1 <- fread(paste0(enrichment_dir,files_interest[f1]) )
      dat1 <- dat1[dat1$gene %in% group_keep, ]
      d1e <-  scale(dat1$enrichment, center = T, scale = T)
      dat2 <- fread(paste0(enrichment_dir, files_interest[f2]) )
      dat2 <- dat2[dat2$gene %in% group_keep, ]
      d2e <- scale(dat2$enrichment, center =T, scale =T) 
      r = cor(d1e, d2e, method = 'spearman')
      output[f1,f2] <- r   
    }
  }
  output  <- data.frame(output) 
  colnames(output) <- toupper(fixNames(signal_names) )
  rownames(output) <- toupper(fixNames(signal_names) )
  output <- cleanDF(output) #fixes names
  output_final[[k]] <- output
}

# output final files
write.table(output_final[[1]], paste0(output_dir, 'arid1a_factors_cor.tab'), col.names = T, row.names = T, sep = '\t') 
write.table(output_final[[2]], paste0(output_dir, 'arid1b_factors_cor.tab'), col.names = T, row.names = T, sep = '\t') 
write.table(output_final[[3]], paste0(output_dir, 'arid2_factors_cor.tab'), col.names = T, row.names = T, sep = '\t') 
write.table(output_final[[4]], paste0(output_dir, 'multi_factors_cor.tab'), col.names = T, row.names =T, sep = '\t')

# those files are useful for clustering the similarity between peaks and assigning the different
# histones/txnfactors to classes
groupnames <- c('Arid1a', 'Arid1b', 'Arid2', 'Multi')

for (k in 1:length(groups)) { 
  out <- NULL
  g <- groups[k]
  files_interest <- getGroupFiles(g, files) 
  files_interest <- files_interest[!files_interest %in% dropname]
  signal_names <- sapply(files_interest, getSigname) 
  for (f1 in 1:length(files_interest)){ 
    dat <- fread(paste0(enrichment_dir, files_interest[f1])) 
    setnames(dat, c('gene', paste(sample(letters, 26, replace =T), collapse='') ) )
    if (!is.null(out)) { 
      out <- out %>% left_join(dat, by = 'gene')
    }
    else{ 
      out <- dat
    }
  
  }
  setnames(out, c('gene', toupper(fixNames(signal_names) )) )
  
  out <- out[,!grepl(toupper('Pol2Forskln*|Srebp1Insln*|Pgc1aForskln*|Mafk_sc*|Hsf1Forskln*|Grp20Forskln*|ErraForskln*|RxlchPcr1x*|RxlchPcr2x*|RxlchV0416101*|RxlchV0422111*'), colnames(out) ), with =F]
  write.table(out, paste0(output_dir, groupnames[k], '_enrichments_bypeaks.csv', sep =''), col.names =T, row.names = F, sep = ',')
}

