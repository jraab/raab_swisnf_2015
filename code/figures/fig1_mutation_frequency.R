# Get mutation data for an arbitrary set of genes. 
# Output will be a plot with the cumulative frequency of mutations
library(cgdsr)
library(ggplot2)
library(dplyr)
library(magrittr)
library(tidyr) 
library(parallel) 
library(RColorBrewer)
library(gridExtra) 
source('code/util/theme_paper.R')
cols <- c('#ef3b2c', '#67a9cf', '#7fbf7b',  '#737373') 
# Set up cgds object and get information about cancers
mycgds <- CGDS('http://www.cbioportal.org/public-portal/') 
all_studies <- getCancerStudies(mycgds) 
genes <- c('ARID1A', 'ARID1B', 'ARID2') %>% toupper() # just to clean up

# Function Defs
get_mut_prof_id <- function(study_id, cgds_obj) { 
  # return profile name for a cancer
  pi <- getGeneticProfiles(cgds_obj, study_id) %>% 
    select(genetic_profile_id) %>% unlist(as.list(.) ) 
  return(pi[grepl(paste(study_id, 'mutation*', sep='_'), pi, perl = T)]) 
}

get_caselist_id <- function(study_id, cgds_obj) { 
  cl <- getCaseLists(cgds_obj, study_id) %>% 
    select(case_list_id) %>% unlist(as.list(.))
  return(cl[grepl(paste(study_id, 'sequenced', sep = '_'), cl, perl = T)] )
}

grab_mutations <- function(study_id, cgds_obj, genes) { 
  # Use a df of studies to return mutation data for each study
  caselist <- get_caselist_id(study_id, cgds_obj) 
  profid <- get_mut_prof_id(study_id, cgds_obj) 
  mutations <- getMutationData(cgds_obj, caselist, profid, genes) 
  return(mutations) 
}

get_casecounts <- function(prof_id, cgds_obj) { 
  study_id <- gsub('_mutation[s]', '', prof_id) 
  cases <- getCaseLists(cgds_obj, study_id) 
  print(cases)
  cases <- cases[grepl(paste(study_id, 'sequenced', sep='_'), cases$case_list_id, perl=T), ]
  return(length(unlist(strsplit(cases$case_ids, '\\s') ) ) )
  }


tally_mutations <- function(mut_df, cgds_obj) {
  mut_counts <- mut_df %>% 
    select(gene_symbol, case_id, genetic_profile_id, mutation_type) %>% 
    filter(!grepl('^Error.*', genetic_profile_id, perl =T) ) %>% 
    select(-genetic_profile_id) %>%
    group_by(case_id, genetic_profile_id) %>% 
    summarise(multiple = ifelse(length(unique(gene_symbol)) > 1, 1, 0), 
              names = gene_symbol) %>% 
    mutate(mutations = ifelse(multiple == 1, 'Multiple', names) ) %>% 
    group_by(mutations, genetic_profile_id) %>% 
    summarise(count = n() ) 
} 

##############################################################################
# Get mutation info
x <- mclapply(all_studies[,1], function(x) try(grab_mutations(x, mycgds, genes) ), mc.cores = 2)  
xx <- do.call('rbind', x) 

mutation_counts <- xx %>% 
select(gene_symbol, case_id, genetic_profile_id, mutation_type) %>% 
  filter(!grepl('^Error.*', genetic_profile_id, perl =T) ) %>% 
  group_by(case_id, genetic_profile_id) %>% 
  summarise(multiple = ifelse(length(unique(gene_symbol)) > 1, 1, 0), 
            names = gene_symbol) %>%
  mutate(mutations = ifelse(multiple == 1, 'Multiple', names) ) %>% 
  group_by(mutations, genetic_profile_id) %>% 
  summarise(count = n() ) 

mutation_counts$casenums <- unlist(mclapply(mutation_counts$genetic_profile_id, function(x) get_casecounts(x, mycgds), mc.cores = 2 ))

pal <- c('#ef3b2c', '#4292c6', '#9ecae1', '#737373')
pal2 <- c(rev(brewer.pal(3, 'Set2')), '#737373') 
mutdf <- mutation_counts %>% 
  mutate(percent = count/casenums * 100) %>%
  group_by(genetic_profile_id) %>% 
  mutate(total = sum(percent) ) %>% 
  ungroup() %>%
  arrange(desc(total)) %>% 
  filter(total > 5) %>%
  mutate(genetic_profile_id = factor(gsub('_mutations', '', genetic_profile_id),
                                     levels = unique(gsub('_mutations', '', 
                                                          genetic_profile_id) ) ) ) 

plt <- mutdf %>% 
  ggplot(aes(x = genetic_profile_id, y = percent, fill = mutations)) + geom_bar(position = 'stack', stat = 'identity') +
    geom_bar(position = 'stack', stat = 'identity', color = 'grey30', show_guide = FALSE) + 
    theme_paper() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5 ), 
                          panel.grid.major.x = element_blank(), 
                          panel.grid.minor.y = element_line(linetype = 'dashed'), 
                          axis.ticks.x = element_blank()) + 
  scale_y_continuous(expand = c(0,0) ) +                
  scale_fill_manual(values = cols) + xlab('') + ylab('Percent Mutated') 

plt

output_df <- mutdf %>% 
  left_join(all_studies, by = c('genetic_profile_id' = 'cancer_study_id') )
 
write.table(output_df, 'output/tables/mutation_data.csv', sep = ',', row.names =F, col.names =T)
ggsave('figures/fig1_mutation_frequency.pdf', plt)
