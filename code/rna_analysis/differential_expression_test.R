# Differential expression testing of swi/snf RNAseq experiments. 
library(DESeq2) 
library(ggplot2) 
library(RColorBrewer)
library(gplots)
library(tidyr)
library(dplyr)
# ------------------------------------------------------------------
# Directories and files that are hard-coded
#proj <- '/magnuson-lab/jraab/analysis/swi_snf_final' # on codon
#proj <- '~/proj/swi_snf_final/' # on linux/mac
counts <- read.csv('data/rna/processed/counts/htseq_counts_output_merged.csv') 
#remove the last 5 rows which contain htseq-coutn info
counts <- counts[1:(nrow(counts)-5), ]
table_output <- 'output/diffexp/tables/' 
plot_output <-  'output/diffexp/plots/' 
if (!file.exists(table_output) & !file.exists(plot_output) ) { 
  dir.create(plot_output, recursive = T) 
  dir.create(table_output, recursive = T) 
}
gene_lengths <- read.table('~/annotations/gencode.v16.gene_lengths')
# Set up experiment design 
# --------------------------------------------------------------------------------
condition <- factor(c('Arid2', 'Arid1a', 'NTG', 'Arid1a', 'Arid2', 'Arid1b', 'Arid2', 'Arid2', 'NTG', 
                    'NTG', 'Arid2', 'Arid1b', 'Arid1a', 'Arid1a', 'NTG', 'NTG', 'Arid1a', 'Arid1b', 'NTG'), 
                    levels = c('NTG', 'Arid1a', 'Arid1b', 'Arid2') )
lane <- factor(c(2,2,3,1,2,
                 3,2,1,1,2,
                 1,3,2,2,1,
                 2,1,3,3) )

genes <- counts[,1]
counts <- counts[,2:ncol(counts)]
row.names(counts) <- genes
coldata <- data.frame(sample = colnames(counts), condition, lane ) 
cds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition + lane )
colnames(cds) <- colnames(counts)
cds<- estimateSizeFactors(cds) 

# Filter genes with low overall counts (Independent Filter Step) 
# -----------------------------------------------------------------------------
min_row_mean <- 2
use <- rowMeans(counts(cds, norm = T) ) > 2 
print(table(use) ) 
cdsFilt <- cds[use, ]

# Determine dispersion and do qc plotting
# -----------------------------------------------------------------------------
cdsFilt <- estimateDispersions(cdsFilt) 
png(paste0(plot_output, 'qc_dispersion_estimate.png') ) 
plotDispEsts(cdsFilt) 
dev.off()


# Differential expression testing
# -----------------------------------------------------------------------------
res <- nbinomWaldTest(cdsFilt) 
arid1a.res <- as.data.frame(results(res, contrast = c('condition', 'Arid1a', 'NTG'), cooksCutoff =F, pAdjustMethod = 'BH' ) )
arid1b.res <- as.data.frame(results(res, contrast = c('condition', 'Arid1b', 'NTG'), cooksCutoff =F , pAdjustMethod ='BH' ))
arid2.res  <- as.data.frame(results(res, contrast = c('condition', 'Arid2' , 'NTG'), cooksCutoff =F , pAdjustMethod ='BH' ) )

# collect total numbers of differentail expressed genes for qc
# -------------------------------------------------------------------------------
diff_gene_totals <- data.frame(gene = c('Arid1a', 'Arid1b', 'Arid2'), 
                               number = c(nrow(na.omit(arid1a.res[arid1a.res$padj < 0.1, ]) ), 
                                          nrow(na.omit(arid1b.res[arid1b.res$padj < 0.05, ]) ), 
                                          nrow(na.omit(arid2.res[arid2.res$padj  < 0.05, ] ) ))
                               )
print(diff_gene_totals)
# Note: I'm using a less stringent cutoff for ARID1A, for some reason only ~ 10% of pvalues < 0.05  
# pass at an FDR < 0.05; Arid1b and arid2 that number is ~50%. Therefore I'm relaxing the FDR to 0.1 
# For arid1a. While I don't recover as many genes as the other two, I do pick up an extra ~ 100
# Significant genes, of which I expect 5 to be false, this is an acceptable tradeoff. 

# QC plots
# --------------------------------------------------------------------------------
rld <- rlogTransformation(cdsFilt, blind = T) 
vsd <- varianceStabilizingTransformation(cdsFilt, blind =T) 
p <- plotPCA(rld, intgroup = c('condition') )
ggsave(filename = paste0(plot_output, 'qc_pca.png'), p )

# Export data for future analyses
# -------------------------------------------------------------------------------
countdata <- counts(cds, norm =T)
fake_rpkm <- function(counts, lengths) {
  #order of lengths must be the same as the order of row.names(counts)
  return(counts*1000/lengths)
}
gene_lengths_reorder <- gene_lengths[gene_lengths$gene_name %in% row.names(countdata), ]
# this asserts that the order of both data frames is the same
table(row.names(countdata) == gene_lengths_reorder$gene_name)
fake_rpkm <- as.data.frame((countdata*1000)/gene_lengths_reorder$exon_length)
fake_rpkm$gene <- row.names(fake_rpkm)
fake_rpkm2 <- fake_rpkm %>% gather(Sample, Value,  -gene) %>% separate(col = Sample, into=c('Sample', 'Rep'), sep='_')

# normalized_count_data (tidy format)
write.table(fake_rpkm2, paste0(table_output, 'count_data_rpkm.csv'), sep=',', col.names =T, row.names =F)

# write full differential genes data

write.table(arid1a.res, paste0(table_output, 'arid1a_full_results.csv'), sep = ',', col.names =T, row.names = T)
write.table(arid1b.res, paste0(table_output, 'arid1b_full_results.csv'), sep = ',', col.names =T, row.names = T) 
write.table(arid2.res,  paste0(table_output, 'arid2_full_results.csv'),  sep = ',', col.names =T, row.names = T)


#write only signficant genes 
arid1a.sig <- na.omit(arid1a.res[arid1a.res$padj < 0.1, ] )
arid1b.sig <- na.omit(arid1b.res[arid1b.res$padj < 0.05, ] )
arid2.sig <- na.omit(arid2.res[arid2.res$padj < 0.05, ])

write.table(arid1a.sig, paste0(table_output, 'arid1a_sig_results.csv'), sep=',', col.names =T, row.names = T) 
write.table(arid1b.sig, paste0(table_output, 'arid1b_sig_results.csv'), sep=',', col.names =T, row.names = T) 
write.table(arid2.sig, paste0(table_output, 'arid2_sig_results.csv'), sep=',', col.names =T, row.names = T)
