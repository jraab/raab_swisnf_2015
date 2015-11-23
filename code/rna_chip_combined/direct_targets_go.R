# Go analysis of direct targets
library(ggplot2)
library(doMC)
library(data.table)
library(plyr)
library(dplyr)
library(tidyr)
cores <- 3 #number of cores, change if less
registerDoMC(cores)
parallel=T #easier to deal with this option here

outdir_figs = 'output/diffexp/plots/'
outrid_tabs = 'output/diffexp/tables/'



#########
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

# Import msigdb information from annotations folder
# -----------------------------------------------------------------------------------------------------------
#I stripped these down a bit - look at only the all curated list and canonical pathways 
#together these two categories cover most curated and go-like lists. 
#these were downloaded 5-5-2014

#all curated
c2.all.msig <- readLines('~/annotations/c2.all.v4.0.symbols.gmt')
c2.all.msig <- strsplit(c2.all.msig, '\t')
names(c2.all.msig) <- sapply(c2.all.msig, function(x) x[1]) 
c2.all.msig <- sapply(c2.all.msig, function(x) x[3:length(x) ])

# canonoical pathways
c2.cp.msig <- readLines('~/annotations/c2.cp.v4.0.symbols.gmt')
c2.cp.msig <- strsplit(c2.cp.msig, '\t')
names(c2.cp.msig) <- sapply(c2.cp.msig, function(x) x[1]) 
c2.cp.msig <- sapply(c2.cp.msig, function(x) x[3:length(x) ])

c5.all.msig <- readLines('~/annotations/c5.all.v4.0.symbols.gmt')
c5.all.msig <- strsplit(c5.all.msig, '\t') 
names(c5.all.msig) <- sapply(c5.all.msig, function(x) x[1]) 
c5.all.msig <- sapply(c5.all.msig, function(x) x[3:length(x) ])

all.paths <- c(c2.all.msig, c5.all.msig) # combine the curated data (c2) and the GO terms (c5) 

# k
# Define Arid lists
#-------------------------------------------------------------------
uni <- 45956 #using this number from the msigdb as universe
direct <- read.csv('output/diffexp/tables/direct_targets.csv')

arid1a <- direct %>% filter(arid == 'Arid1a') %>% select(rowname) %>% as.list(.)
arid1b <- direct %>% filter(arid == 'Arid1b') %>% select(rowname) %>% as.list(.) 
arid2  <- direct %>% filter(arid == 'Arid2')  %>% select(rowname) %>% as.list(.)

allp_hyperg <- function(x, ...){ 
  p <- get.hyperg(x, all.paths, uni, ntop =30, parallel =T, ...)
  return(p)
  }
go_hyperg <- function(x, ...) { 
  p <- get.hyperg(x, c5.all.msig, uni, ntop =30, parallel =T, ...) 
  return(p)
  }

arid1a_allp <- allp_hyperg(arid1a$rowname, returnplot =T)
go_hyperg(arid1a$rowname)
ggsave('output/diffexp/plots/arid1a_allmsigdirect.pdf', arid1a_allp)
arid1b_allp <- allp_hyperg(arid1b$rowname, returnplot =T)
go_hyperg(arid1b$rowname) 
ggsave('output/diffexp/plots/arid1b_allmsigdirect.pdf', arid1b_allp)
arid2_allp <- allp_hyperg(arid2$rowname, returnplot = T)
go_hyperg(arid2$rowname)
ggsave('output/diffexp/plots/arid2_allmsigdirect.pdf', arid2_allp)

# intersect the arid1b and arid2 lists - and separate up vs down regulated
arid1b_arid2 <- intersect(arid1b$rowname, arid2$rowname) 
allp_hyperg(arid1b_arid2)
go_hyperg(arid1b_arid2) 

# get log2fold change for thse 
arid1b_arid2_direct <- direct %>% 
  filter(rowname %in% arid1b_arid2) %>%
  filter(!arid == 'Arid1a') %>%
  spread(arid, log2FoldChange) 

arid1b_arid2_direct %>%
  ggplot(aes(x=Arid1b, y =Arid2)) + 
  geom_point(size =2) + 
  theme_classic() + 
  geom_hline(x=0, color = 'grey30', linetype='dashed') + 
  geom_vline(y=0, color = 'grey30', linetype='dashed') + 
  theme(axis.text = element_text(size = 20), 
        axis.title = element_text(size =24) ) + 
  xlab('Log2FC siArid1b/siNTG') + ylab('Log2FC siArid2/siNTG') -> arid1b_arid2_directfc
arid1b_arid2_directfc  
ggsave('output/diffexp/plots/arid1b_arid2_direct_log2fc.pdf', arid1b_arid2_directfc)


arid1b_arid2_repress <- arid1b_arid2_direct %>% 
  filter(Arid1b > 0 & Arid2 > 0 ) 

arid1b_arid2_activate <- arid1b_arid2_direct %>% 
  filter(Arid1b < 0 & Arid2 < 0 )

#compare go and msigdb
arid1b_arid2_rep_msig_p <- allp_hyperg(arid1b_arid2_repress$rowname, returnplot=T)
go_hyperg(arid1b_arid2_repress$rowname)
ggsave('output/diffexp/plots/arid1b_arid2_direct_repressed_msigcat.pdf', arid1b_arid2_rep_msig_p)

arid1b_arid2_act_msig_p <- allp_hyperg(arid1b_arid2_activate$rowname, returnplot =T) 
go_hyperg(arid1b_arid2_activate$rowname) 
ggsave('output/diffexp/plots/arid1b_arid2_direct_activated_msigcat.pdf', arid1b_arid2_act_msig_p)

#create a bed file that is just the arid1b_arid2_repressed_genes
all_genes <- read.csv('~/annotations/gencode.v16.annotation.bed', header=F, sep='\t')
output <- all_genes[all_genes$V4 %in% arid1b_arid2_repress$rowname,]
write.table(output, 'output/diffexp/tables/arid1b_arid2_direct_repressed.bed', row.names= F, col.names =F, quote =F, sep ='\t')


