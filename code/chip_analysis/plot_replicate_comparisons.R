setwd('/magnuson-lab/jraab/analysis/swi_snf_final/')
inputdir <- 'data/chip/processed/replicate_coverages/'
outdir <- 'output/'
library(ggplot2) 
#takes a stub name name for chipseq and makes a plot of the correlation between the two 
args <- commandArgs(TRUE)
print(args)
file1 <- read.table( paste0(inputdir, args[1], '_rep1.10kbcov.txt'), header = F , sep='\t') 
file2 <- read.table( paste0(inputdir, args[1], '_rep2.10kbcov.txt'), header = F, sep='\t') 

print ('files in')
df <- data.frame( Rep1 = file1$V4, Rep2 = file2$V4) 
print (head(df) )
#normalize by total count in each column which should e a good proxy for mapped reads
df$Rep1 <- log2((df$Rep1+1)/sum(df$Rep1+1) ) * 1e6 # reads per million
df$Rep2 <- log2((df$Rep2+1)/sum(df$Rep2+1) ) * 1e6  
corstat <- cor(df$Rep1, df$Rep2)
g <- ggplot(df, aes(x = Rep1, y = Rep2 ) )+ geom_point() + geom_text(data=data.frame(x = max(df$Rep1) * 0.9, y = min(df$Rep2), label=corstat), aes(label=label, x=x, y=y) )    
ggsave(paste0(outdir,args[1], '10kbcv_cor.pdf'),  g) 
                  


