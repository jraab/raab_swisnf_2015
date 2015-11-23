# Compare TPM data for ENCODE cell lines to answer reviewer question about 
# ARID expression levels in my cell type compared to others. 
library(readr) 
library(dplyr) 
library(tidyr) 
library(biomaRt) 
library(ggplot2) 
source('code/util/theme_paper.R') 

# Downloaded all data from ENCODE (2015/10/11)
files <- list.files('data/external/encode_expr', '^ENC.*.tsv', full.names = T, )
metadata <- 'data/external/encode_expr/metadata.tsv'
metadata <- read_tsv(metadata) 
gene_data <- metadata %>% filter(`Output type` == 'gene quantifications') 
gene_data_acc <- gene_data$`File accession`

#creat a mapping of cell type and accession number 
gene_data_min <- gene_data %>% dplyr::select(`File accession`, `Biosample term name`) 

# Get expression data for the gene of interest
goi <- c('ARID1A', 'ARID1B', 'ARID2') 
# Map these to ensembl IDs
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'www.ensembl.org')
mart_res <- getBM(attributes = c('ensembl_transcript_id', 'ensembl_gene_id', 'hgnc_id', 'external_gene_name'), mart = mart) 
gene_mapping <- mart_res %>% 
  filter(external_gene_name %in% goi) 

# Create a data frame of values for just the correct genes
get_expr_info <- function(genelist, expr_files){ 
  l <- list()
  for (i in 1:length(expr_files)) { 
    fname <- basename(expr_files[i])
    fname <- gsub('.tsv', '', fname)
    d <- read_tsv(expr_files[i], col_types = 'ccddddddddddddd') 
    d <- d %>% separate(gene_id, into = c('gene_id', 'leftover.dot'), sep = '\\.', extra = 'drop') 
    d <- d %>% filter(gene_id %in% genelist) 
    l[[i]] <- data.frame(gene_id = d$gene_id, TPM= d$TPM, FPKM = d$FPKM, exp_name = fname)
    l[[i]]$exp_name <- fname
  }
  out <- do.call(rbind.data.frame, l)
  
  return(out) 
  }

# # need to match file names to gene quantifications
# gene_data_files <- list.files('data/external/encode_expr', pattern = paste(gene_data_acc, sep = '|', collapse='|'), full.names = T ) 
# test <- get_expr_info(unique(gene_mapping$ensembl_gene_id), gene_data_files) 
# 
# # Now merge the experiment names back into the cell type names
# # 
# g <- gene_mapping %>% 
#   dplyr::select(-ensembl_transcript_id, -hgnc_id) %>%
#   unique(.) %>%
#   left_join(test, by = c('ensembl_gene_id'='gene_id')) %>% 
#   left_join(gene_data_min, by = c('exp_name'= 'File accession') ) %>% 
#   group_by(`Biosample term name`, external_gene_name) %>%
#   summarise(mean_tpm = mean(TPM, na.rm = T) ) %>%
#   ungroup() %>%
#   arrange(desc(mean_tpm)) %>%
#   mutate(bsample =factor(`Biosample term name`, levels = unique(as.character(`Biosample term name`) ))) 
# ggplot(g, aes(x= bsample, y = mean_tpm)) + geom_point() + facet_wrap(~external_gene_name) + 
#   geom_point(data = g[g$bsample == 'HepG2',], color = 'red2', size = 3) + coord_flip()
# 
# # All data points
# g2 <- gene_mapping %>% 
#   dplyr::select(-ensembl_transcript_id, -hgnc_id) %>%
#   unique(.) %>%
#   left_join(test, by = c('ensembl_gene_id'='gene_id')) %>% 
#   left_join(gene_data_min, by = c('exp_name'= 'File accession') ) 
#  
# # Data frame summarized by cell line
# g3 <- g2 %>% 
#   group_by(`Biosample term name`, external_gene_name) %>% 
#   summarise(mean_tpm = mean(TPM, na.rm = T) ) %>%
#   ungroup() %>%
#   arrange(mean_tpm) %>% 
#   mutate(`Biosample term name` = factor(`Biosample term name`, 
#                                         levels = as.character(unique(`Biosample term name`))))
# 
# # Order the full data frame by the levels in g3
# g2 <- g2 %>% 
#   mutate(`Biosample term name` = factor(`Biosample term name`, levels = levels(g3$`Biosample term name`)))
#   
# plt <- ggplot(g2, aes(x = `Biosample term name`, y= TPM)) + geom_point() + facet_wrap(~external_gene_name) + 
#     geom_point(data = g3, aes(x = `Biosample term name`, y = mean_tpm), shape =108 , color = 'red2', size = 5) + coord_flip() + theme_paper() + xlab('') + theme(panel.border = element_rect(color = 'grey30', fill = NA), 
#                                                                                                                                                                  axis.line = element_blank(), panel.grid.major = element_line(color = 'grey90'))
# 
# ggsave('figures/supp_tpm_bycellline.pdf', plt) 

# Could alternatively filter to use only certain files - 
gene_data <- metadata %>% filter(`Output type` == 'gene quantifications', `Biosample subcellular fraction term name` == '')  
gene_data_acc <- gene_data$`File accession`
gene_data_files <- list.files('data/external/encode_expr', pattern = paste(gene_data_acc, sep = '|', collapse='|'), full.names = T ) 
test <- get_expr_info(unique(gene_mapping$ensembl_gene_id), gene_data_files) 
# 
g <- gene_mapping %>% 
  dplyr::select(-ensembl_transcript_id, -hgnc_id) %>%
  unique(.) %>%
  left_join(test, by = c('ensembl_gene_id'='gene_id')) %>% 
  left_join(gene_data_min, by = c('exp_name'= 'File accession') ) %>% 
  group_by(`Biosample term name`, external_gene_name) %>%
  summarise(mean_tpm = mean(TPM, na.rm = T) ) %>%
  ungroup() %>%
  arrange(desc(mean_tpm)) %>%
  mutate(bsample =factor(`Biosample term name`, levels = unique(as.character(`Biosample term name`) ))) 
ggplot(g, aes(x= bsample, y = mean_tpm)) + geom_point() + facet_wrap(~external_gene_name) + 
  geom_point(data = g[g$bsample == 'HepG2',], color = 'red2', size = 3) + coord_flip()

# All data points
g2 <- gene_mapping %>% 
  dplyr::select(-ensembl_transcript_id, -hgnc_id) %>%
  unique(.) %>%
  left_join(test, by = c('ensembl_gene_id'='gene_id')) %>% 
  left_join(gene_data_min, by = c('exp_name'= 'File accession') ) 

# Data frame summarized by cell line
g3 <- g2 %>% 
  group_by(`Biosample term name`, external_gene_name) %>% 
  summarise(mean_tpm = mean(TPM, na.rm = T) ) %>%
  ungroup() %>%
  arrange(mean_tpm) %>% 
  mutate(`Biosample term name` = factor(`Biosample term name`, 
                                        levels = as.character(unique(`Biosample term name`))))

# Order the full data frame by the levels in g3
g2 <- g2 %>% 
  mutate(`Biosample term name` = factor(`Biosample term name`, levels = levels(g3$`Biosample term name`)))

plt <- ggplot(g2, aes(x = `Biosample term name`, y= TPM)) + geom_point() + facet_wrap(~external_gene_name) + 
  geom_point(data = g3, aes(x = `Biosample term name`, y = mean_tpm), shape =108 , color = 'red2', size = 5) + 
  coord_flip() + theme_paper() + xlab('') + 
  theme(panel.border = element_rect(color = 'grey30', fill = NA),
        axis.line = element_blank(), panel.grid.major = element_line(color = 'grey90'))
ggsave('figures/supp_tpm_bycellline.pdf', plt) 

