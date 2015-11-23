# This script was an attempt to answer a reviewer question about eRNA transcription in our knockdowns# However, we used poly A select mRNA to generate libraries, and it does not appear we observe any 
# Strong evidence of transcription at enhancers. 
# Nevertheless we attempted to quantitate the level of transcription at 'enhancers' and see if 
# any differences were observed in our various ARID knockdowns. 
# See code/rna_analysis/count_rna_at_enhancers.sh for the htseq-count script
# Below we check for differential expression - but only really see either technical noise 
# (from duplicated reads)
# Or transcription that likely is associated with gene (near exon etc).


library(ggplot2) 
library(tidyr) 
library(dplyr) 
library(DESeq2) 
library(readr) 
source('code/util/chippeak_annotation.R')

active <- read_csv('data/rna/processed/counts/enhancers/htseq_counts_output_merged_active.csv')
active <- active %>% filter(!id %in% c('alignment_not_unique', 'not_aligned', 'too_low_aQual', 'ambiguous') )
inactive <- read_csv('data/rna/processed/counts/enhancers/htseq_counts_output_merged_inactive.csv')
inactive <- inactive %>% filter(!id %in% c('alignment_not_unique', 'not_aligned', 'too_low_aQual', 'ambiguous') )

active_totals <- colSums(active[2:ncol(active)]) 
inactive_totals <- colSums(inactive[2:ncol(inactive)])

# Set up experimental design - treat active and inactive enhancers completely separately. 
# ------------------------------------------------------------------------------------------------

active_conditions <- factor(c('ARID1A', 'ARID2', 'NTG', 'ARID1A', 'ARID2', 'ARID2', 'ARID1A', 
                            'ARID1B', 'NTG', 'ARID1A', 'NTG', 'ARID1B', 'NTG', 'ARID2', 'ARID1A', 'ARID2', 'NTG', 'NTG', 'ARID1B'), 
                            levels = c('NTG', 'ARID1A', 'ARID1B', 'ARID2') )
inactive_conditions <- factor(c('NTG', 'NTG', 'NTG', 'NTG', 'ARID1A', 'ARID2', 'ARID2', 'ARID1B', 
                                'NTG', 'ARID1A', 'ARID1A', 'ARID2', 'ARID2', 'ARID1A', 'ARID1A',
                                'ARID1B', 'ARID2', 'NTG', 'ARID1B'),
                              levels = c('NTG', 'ARID1A', 'ARID1B', 'ARID2') )

inactive_ids <- inactive[,1]
inactive <- inactive[,2:ncol(inactive)]
rownames(inactive) <- inactive_ids$id

active_ids <- active[,1]
active <- active[,2:ncol(active)]
row.names(active) <- active_ids$id

active_coldata <- data.frame(sample = colnames(active), condition = active_conditions) 
active_cds <- DESeqDataSetFromMatrix(countData = active, colData = active_coldata, design = ~ condition)
colnames(active_cds) <- colnames(active) 

inactive_coldata <- data.frame(sample = colnames(inactive), condition = inactive_conditions)
inactive_cds <- DESeqDataSetFromMatrix(countData = inactive, colData = inactive_coldata, design = ~ condition)
colnames(inactive_cds) <- colnames(inactive) 

# Estimate size factors for each separately and compare - if not the same need to do something since 
# these came from the same libraries. I kept the numbers of reads mapping to no features 
# This was done as a measure of the mRNA mapping reads - which I didn't think it appropriate to ignore
# Since it is a large fraction of all reads in these libraries (93-97% of reads do not map to 'enhnacers')

active_cds <- estimateSizeFactors(active_cds) 
inactive_cds <- estimateSizeFactors(inactive_cds) 

# Plot them to check relative ratio
active_df_sizes <- as.data.frame(sizeFactors(active_cds)) %>% add_rownames() 
inactive_df_sizes <- as.data.frame(sizeFactors(inactive_cds)) %>% add_rownames() 
df <- left_join(active_df_sizes, inactive_df_sizes, by = 'rowname') 
colnames(df) <- c('rowname', 'active', 'inactive')
df %>% gather(group, val, -rowname) %>% 
  ggplot(aes(x=rowname, y = val, fill = group)) + geom_bar(position = 'dodge', stat = 'identity') 
# Everything looks fine - since the no feature reads dominate - the size factors are pretty close. 
# Can proceed. 

# -------------------------------------------------------------------------------------------------
min_row_mean <-.1
use_act <- rowMeans(counts(active_cds, norm = T)) > min_row_mean
use_in <- rowMeans(counts(inactive_cds, norm = T)) > min_row_mean
print(table(use_act) )
print(table(use_in) )

# Filters a lot out (7431 active (37%), 697 inactive (6.5%)) - ratio makes sense here, since inactive should not 
# making eRNA.
act_filt <- active_cds[use_act,]
in_filt  <- inactive_cds[use_in, ]

# Determine dispersion and do qc plotting
# -----------------------------------------------------------------------------
act_filt <- estimateDispersions(act_filt) 
in_filt  <- estimateDispersions(in_filt) 

#png(paste0(plot_output, 'qc_dispersion_estimate.png') ) 
plotDispEsts(act_filt)  
plotDispEsts(in_filt) 

# DE calling
# --------------------------------------------------------------------------------------------------
act_test <- nbinomWaldTest(act_filt)  

act_res_arid1a <- as.data.frame(results(act_test, contrast = c('condition', 'ARID1A', 'NTG'), cooksCutoff =F, pAdjustMethod = 'BH' ) ) %>%
  add_rownames() %>%
  mutate(ARID='ARID1A', group = 'Active')
act_res_arid1b <- as.data.frame(results(act_test, contrast = c('condition', 'ARID1B', 'NTG'), cooksCutoff =F, pAdjustMethod = 'BH' ) ) %>%
  add_rownames() %>%
  mutate(ARID='ARID1B', group = 'Active')
act_res_arid2 <- as.data.frame(results(act_test, contrast = c('condition', 'ARID2', 'NTG'), cooksCutoff =F, pAdjustMethod = 'BH' ) )  %>%
  add_rownames() %>%
  mutate(ARID='ARID2', group = 'Active') 

inact_test <- nbinomWaldTest(in_filt) 
inact_res_arid1a <- as.data.frame(results(inact_test, contrast = c('condition', 'ARID1A', 'NTG'), cooksCutoff =F, pAdjustMethod = 'BH' ) ) %>%
  add_rownames() %>%
  mutate(ARID = 'ARID1A', group = 'Inactive') 
inact_res_arid1b <- as.data.frame(results(inact_test, contrast = c('condition', 'ARID1B', 'NTG'), cooksCutoff =F, pAdjustMethod = 'BH') ) %>%
  add_rownames() %>%
  mutate(ARID = 'ARID1B', group = 'Inactive') 
inact_res_arid2 <- as.data.frame(results(inact_test, contrast = c('condition', 'ARID2', 'NTG'), cooksCutoff =F, pAdjustMethod = 'BH') )  %>%
  add_rownames() %>%
  mutate(ARID = 'ARID2', group = 'Inactive')
                            
# Simple table of numbers of changed genes
# --------------------------------------------------------------------------------------------------

all_data <- rbind(act_res_arid1a, act_res_arid1b, act_res_arid2, 
                  inact_res_arid1a, inact_res_arid1b, inact_res_arid2) 

all_data %>%
  filter(padj < 0.1) %>%
  mutate(direction = ifelse(log2FoldChange < 0, 'Down', 'Up')) %>%
  group_by(ARID, group, direction) %>% 
  summarise(count = n(), ) %>%
  ggplot(aes(x= ARID, y = count, fill = direction)) + geom_bar(position = 'dodge', stat = 'identity') + 
  facet_wrap(~group)

# Checking manually by hand - I see most places as exons that are bleeding into enahncer (or vice versa) .
# These are unlikely to be real enhancer transcription - to demonstrate I'll filter my list against exons from gencode. 

gencode <- gtf2gr('data/external/gencode.v16.annotation.gtf')
combined_data <- all_data %>% 
  filter(!rowname == 'no_feature') %>%
  separate(rowname, into = c('start', 'end'), sep = '-', remove = F) %>%
  separate(start, into = c('chr', 'start'), sep = ':')

combined_gr <- with(combined_data, GRanges(chr, IRanges(as.numeric(start),as.numeric(end)), name = rowname))

# Filter out any of these enhancers that overlap an exon. 
filtered <- combined_gr[!overlapsAny(combined_gr, gencode[gencode$type == 'exon',])]
# This drops about 6000 'enhancers'
# Merge this information back in 

combined_data$inExon <- ifelse(combined_data$rowname %in% filtered$name, 0, 1) 
combined_data$isSig <- ifelse(combined_data$padj < 0.1, 1, 0) 

# Now tally whether my 'signficantly' changed 'enhancers' are just exons from altered genes
all_data_byexon <- combined_data %>% 
  filter(padj < 0.1) %>%
  group_by(ARID, group) %>%
  summarise(percent = sum(inExon)/n(), count_in_exons = sum(inExon)) 

# We find a total of 137 'enhancers' that change' - but the majority of these
# are found in exons
sum(all_data_byexon$count_in_exons) 

# We can look at those at those that are not 
combined_data %>% 
  filter(padj < 0.1, inExon == 0) %>%
  View

#This give us 78 'enhancers' - there is not strong evidence of transcription here, and mostly likely 
# this is just noise from small levels of duplication in the library. 