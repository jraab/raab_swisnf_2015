#!/usr/bin/env/python
################
#general script to take a bam file and make a bigwig file
# do not store intermediates - since these take pu stupid amounts of space
#program requires the ucsc utility bedgraphToBigWig to be somewhere in the $PATH
#Usage:bamToWig.py -b INPUT.bam -g hg19|mm9 -s FALSE --norm  
# make sure module load samtools is in submit script
###################
import argparse
import HTSeq
import numpy
import os
import sys

parser = argparse.ArgumentParser(description = 'Convert BAM files to wiggle tracks' ) 
parser.add_argument('-b', help='input bam' )
parser.add_argument('--norm', help = 'Normalize by read depth' , action='store_true')
parser.add_argument('-s', help = 'separate files for + and - strand', action='store_true', default=False)
parser.add_argument('-g', help = 'genome name e.g hg19')
parser.add_argument('-d', help = 'outputdirectory') 
parser.add_argument('-frag', help = 'fragment size of library', default = 200) 
args = parser.parse_args()

binsize = 30 # looksl like encode bins things at about 30 bases or so)
home = 'magnuson-lab/jraab/'
################################
def calc_coverage(bamfile, chrom_dict, fragmentsize, stranded=False):
	read_count =0 
	millions =1
	bamfile = HTSeq.BAM_Reader(bamfile)
	fragmentsize = fragmentsize #smooth out reads based on fragmentsizes	
	if stranded==True: 
		cvg = HTSeq.GenomicArray(chrom_dict, stranded=True, typecode='d', storage='step')
	else: 
		cvg = HTSeq.GenomicArray(chrom_dict, stranded=False, typecode='d', storage='step')	
	for almnt in bamfile: 
		if almnt.aligned: 
			end_of_chrom = chrom_dict[almnt.iv.chrom]	
			almnt.iv.length = fragmentsize	
			if almnt.iv.start <1: 
				almnt.iv.start =1
			if almnt.iv.end > end_of_chrom: 
				almnt.iv.end = end_of_chrom
			cvg[almnt.iv]+=1
			read_count += 1 
			if read_count % 1e6 == 0:
				processed = int(1e6*millions)
				print ' Processed %.1E reads' % processed

				millions +=1
	if args.norm: 
		normfactor = 1e6*binsize/read_count
		for iv, val in cvg.steps(): 
			cvg[iv].apply(lambda x: x* normfactor)	
	else: 
		cvg = cvg 
	return(cvg)	

def read_genome(genome): 
	with open(genome) as g:
		chrom_sizes=dict() 
		for l in g:	
			if l.startswith('chrom'): 
				next
			else :
				chrom_sizes[l.split('\t')[0]] = int(l.split('\t')[1])
	return(chrom_sizes) 
###############################
#this file must exist
if args.g == 'hg19': 
   genome = '/magnuson-lab/shared/jraab/genome_sizes/genome_hg19.txt'
   chrom_sizes = read_genome(genome)
elif args.g =='mm9':
   chrom_sizes = '/magnuson-lab/shared/jraab/genome_sizes/genome_mm9.txt'
   chrom_sizes = read_genome(genome)
else: 
   sys.exit('Need to input a valid genome hg19 or mm9')


name = args.b.split('/')[-1]
name = name.split('.')[0]
#change this if you want output to go somewhere else
dir = args.d

if not os.path.exists(dir): 
   os.makedirs(dir)

filename=os.path.basename(args.b).split('.')[0]
if args.s == True:	
	#for RNAseq stranded protocol these will come out backwards. 
	genome_cvg= calc_coverage(args.b, chrom_sizes, stranded=True)
	genome_cov_plus.write_bedgraph(dir+filename+'_plus.bg', strand='+')
	genome_cov_minus.write_bedgraph(dir+filename+'_minus.bg', strand='-')
else: 
	genome_cov= calc_coverage(args.b, chrom_sizes, stranded=False, fragmentsize = int(args.frag)) 
	genome_cov.write_bedgraph_file(dir+filename+'.bg')

