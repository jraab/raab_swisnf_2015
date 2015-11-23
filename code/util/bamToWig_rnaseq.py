#!/usr/bin/env/python
################
#general script to take a bam file and make a bigwig file
# do not store intermediates - since these take pu stupid amounts of space
#program requires the ucsc utility bedgraphToBigWig to be somewhere in the $PATH
import argparse
import HTSeq
import numpy
import pybedtools
import os

parser = argparse.ArgumentParser(description = 'Convert BAM files to wiggle tracks' ) 
parser.add_argument('-b', help='input bam' )
parser.add_argument('--norm', help = 'Normalize by read depth' , action='store_true')
parser.add_argument('-s', help = 'separate files for + and - strand', action='store_true', default=False)
parser.add_argument('-g', help = 'genome name e.g hg19')
parser.add_argument('-d', help = ' output director') 
args = parser.parse_args()

################################
def calc_norm_factor(bamfile):	
	bam_reader = HTSeq.BAM_Reader(args.b)
        read_counter = 0
	for read in bam_reader: 
		if read.aligned: 
			read_counter +=1
	return( 1e6 /read_counter) 

###############################
#this file must exist
if args.g == 'hg19': 
   genome = '/magnuson-lab/jraab/annotations/genome_hg19.txt'
elif args.g =='mm9':
   chrom_sizes = '/magnuson-lab/shared/jraab/annotations/genome_mm9.txt'
else: 
   sys.exit('Need to input a valid genome hg19 or mm9')

name = args.b.split('/')[-1]
name = name.split('.')[0]
dir = args.d

if not os.path.exists(dir): 
   os.makedirs(dir)

if args.norm:
	normfactor = calc_norm_factor(args.b)

bt = pybedtools.BedTool(args.b)
#hack below because I'm using strand specific rnaseq (where the read comes from the reverse
#TODO add formal option for specifying strnad and name output correctly. For now all strands below
# are reversed 
if args.s:
        output_plus = name+'_plus.bw'
        output_minus = name+'_minus.bw'
        if args.norm:
		bedgraph_plus = bt.genome_coverage(genome=args.g, split=True, bg=True, strand="-", scale=normfactor)
		bedgraph_minus =bt.genome_coverage(genome=args.g, split=True, bg=True, strand='+', scale=normfactor)
	else:
		bedgraph_plus = bt.genome_coverage(genome=args.g, split=True, bg=True, strand="-")
		bedgraph_minus = bt.genome_coverage(genome=args.g, split=True, bg=True, strand="+")
        cmds = ['bedGraphToBigWig', bedgraph_plus.fn, genome, dir+output_plus] 
        os.system(' '.join(cmds))
        cmds = ['bedGraphToBigWig', bedgraph_minus.fn, genome, dir+output_minus] 
        os.system(' '.join(cmds))
else:
	output = name+'.bw'
	if args.norm: 
		bedgraph = bt.genome_coverage(genome=args.g, split=True, bg=True, scale=normfactor)
	else:
		bedgraph = bt.genome_coverage(genome=args.g, split=True,  bg=True, scale=normfactor)

	cmds = ['bedGraphToBigWig', bedgraph.fn, genome, dir+output]
        os.system(' '.join(cmds))


