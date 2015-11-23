# GO analysis of ALONE vs JOINTLY Regulated Genes
# Reviewer requested this analysis be doen. 

library(ggplot2)
library(doMC)
library(data.table)
library(plyr)
library(dplyr)
library(tidyr)
library(readr) 
cores <- 2 #number of cores, change if less
registerDoMC(cores)
parallel=T #easier to deal with this option here

outdir_figs = 'output/diffexp/plots/'
outdir_tabs = 'output/diffexp/tables/go/'

if (!dir.exists(outdir_tabs)){ dir.create(outdir_tabs)}

##################################################################################################3
#function defs 
#function to perform hypergeometric test on two lists of genes
hyperg <- function(my.genes, test.genes, universe_size) { 
  #my.genes = genes on my list
  #test.genes = enrichment set to test against
  #universe_size = total number of genes (according to msigdb 45956 - I'll use this) 
  overlap <- length(intersect(my.genes, test.genes))  
  draws <- length(test.genes) 
  p <- phyper(overlap-1, draws , universe_size-draws, length(my.genes), lower.tail = F )
  
  #phyper(#number overlap, #number in gene set, #number in universe - #number gene set, #number on my list)
  return(data.frame(pval=p, Fraction_Overlap=overlap/draws,  numSet=draws))
}

#convienence function to call hyperg on many at once and store the results) 
get.hyperg <- function(my.genes, test.set, universe_size, ntop=30, returnplot=F, cutoff=0.001, parallel=T, ...) { 
  require(ggplot2) 
  require(data.table) 
  if (parallel==T) { 
    require(doMC)}
  pval <- data.table(ldply(test.set, .fun = hyperg,
                           my.genes=my.genes, universe_size=universe_size, .parallel=parallel))
  print(dim(pval) ) 
  pval[,'qval':= p.adjust(pval, method='BH', n=length(pval)),]
  print(dim(pval) ) 
  sig <- pval[qval < cutoff, ]
  sig <- sig[order(qval)]
  top <- sig[1:ntop]
  top <- top[,.id:=factor(.id, levels=rev(.id)),]
  if (nrow(sig) < 1 ) { 
    print ('No pathways of significance found' )
    return()
  }
  else {
    g <- ggplot(na.omit(top), aes(y=-log10(qval), x=.id)) + geom_point(aes(size = Fraction_Overlap), pch = 19 )
    g <- g + theme_classic() +  ylab('-Log10 qVal') +xlab('') + coord_flip() + expand_limits(y =0 )
    print(g)
    if (returnplot==T) { 
      return(g) 
    }
    else {
      return(sig) 
    }
  }
}

# Test for enriched pathways using all pathways
allp_hyperg <- function(x, ...){ 
  p <- get.hyperg(x, all.paths, uni, ntop =30, parallel =T, ...)
  return(p)
}

# Test for enriched pathways using only GO annotations
go_hyperg <- function(x, ...) {
  p <- get.hyperg(x, c5.all.msig, uni, ntop =30, parallel =T, ...) 
  return(p)
}

# Import msigdb information from annotations folder
# -----------------------------------------------------------------------------------------------------------
#I stripped these down a bit - look at only the all curated list and canonical pathways 
#together these two categories cover most curated and go-like lists. 
#these were downloaded 5-5-2014

#all curated
c2.all.msig <- readLines('data/external/go_anno/c2.all.v4.0.symbols.gmt')
c2.all.msig <- strsplit(c2.all.msig, '\t')
names(c2.all.msig) <- sapply(c2.all.msig, function(x) x[1]) 
c2.all.msig <- sapply(c2.all.msig, function(x) x[3:length(x) ])

# canonoical pathways
c2.cp.msig <- readLines('data/external/go_anno/c2.cp.v4.0.symbols.gmt')
c2.cp.msig <- strsplit(c2.cp.msig, '\t')
names(c2.cp.msig) <- sapply(c2.cp.msig, function(x) x[1]) 
c2.cp.msig <- sapply(c2.cp.msig, function(x) x[3:length(x) ])

c5.all.msig <- readLines('data/external/go_anno/c5.all.v4.0.symbols.gmt')
c5.all.msig <- strsplit(c5.all.msig, '\t') 
names(c5.all.msig) <- sapply(c5.all.msig, function(x) x[1]) 
c5.all.msig <- sapply(c5.all.msig, function(x) x[3:length(x) ])

all.paths <- c(c2.all.msig, c5.all.msig) # combine the curated data (c2) and the GO terms (c5) 

# Define gene lists 
# -------------------------------------------------------------------------------------------------
uni <- 45956 #using this number from the msigdb as universe
#uni <- 17176 #Just the genes we measure expression for - should filter the annotations to reflect this list
# using the larger set, because in this method if these genes were important we could have detected them
# Therefore lack of expression should have meaning here. 

arid1a <- read_csv('output/diffexp/tables/arid1a_diff_genes.csv')
arid1a_alone <- arid1a %>% filter(alone_or_joint == 'Alone') %>% mutate(knockdown = 'ARID1A') 
arid1a_joint <- arid1a %>% filter(alone_or_joint == 'Jointly') %>% mutate(knockdown = 'ARID1A')

arid1b <- read_csv('output/diffexp/tables/arid1b_diff_genes.csv') 
arid1b_alone <- arid1b %>% filter(alone_or_joint == 'Alone') %>% mutate(knockdown = 'ARID1B') 
arid1b_joint <- arid1b %>% filter(alone_or_joint == 'Jointly') %>% mutate(knockdown = 'ARID1B')

arid2 <- read_csv('output/diffexp/tables/arid2_diff_genes.csv')
arid2_alone <- arid2 %>% filter(alone_or_joint == 'Alone') %>% mutate(knockdown = 'ARID2')
arid2_joint <- arid2 %>% filter(alone_or_joint == 'Jointly') %>% mutate(knockdown = 'ARID2')

#Write a supplemental table that has information for all gene lists. 

combined <- rbind(arid1a_alone, arid1a_joint, arid1b_alone, arid1b_joint, arid2_alone, arid2_joint) 
write_csv(combined, 'output/diffexp/tables/all_arids_alone_vs_joint.csv') 

# Test for enrichments
arid1a_alone_allp <- allp_hyperg(arid1a_alone$gene, returnplot =F)
go_hyperg(arid1a_alone$gene)

arid1a_joint_allp <- allp_hyperg(arid1a_joint$gene, retrunplot = F) 
go_hyperg(arid1a_joint$gene)

arid1b_alone_allp <- allp_hyperg(arid1b_alone$gene, returnplot = F) 
go_hyperg(arid1b_alone$gene)

arid1b_joint_allp <- allp_hyperg(arid1b_joint$gene, returnplot = F) 
go_hyperg(arid1b_joint$gene) 

arid2_alone_allp <- allp_hyperg(arid2_alone$gene, returnplot = F) 
go_hyperg(arid2_alone$gene)

arid2_joint_allp <- allp_hyperg(arid2_joint$gene, returnplot = F) 
go_hyperg(arid2_joint$gene) 

write_csv(arid1a_alone_allp, 'output/diffexp/tables/go/arid1a_alone_paths.csv')
write_csv(arid1a_joint_allp, 'output/diffexp/tables/go/arid1a_joint_paths.csv')
write_csv(arid1b_alone_allp, 'output/diffexp/tables/go/arid1b_alone_paths.csv')
write_csv(arid1b_joint_allp, 'output/diffexp/tables/go/arid1b_joint_paths.csv')
write_csv(arid2_alone_allp, 'output/diffexp/tables/go/arid2_alone_paths.csv') 
write_csv(arid2_joint_allp, 'output/diffexp/tables/go/arid2_joint_paths.csv') 

