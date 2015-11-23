#!/usr/bin/env python 
# script to convert merge nearby peaks and filter peaks on the blacklist from macs called chipseq peaks
# will also  fix the silly macs named peaks by lopping off al the directory garbage
###################
import pybedtools as pbt 
import os 
import pandas as pd

directory = '/magnuson-lab/jraab/analysis/swi_snf_final/'
peakdir = directory + 'output/macs_peaks/'
patt = 'broadPeak'
blacklist = '/magnuson-lab/jraab/annotations/combinedblacklist.expanded.bed'
blacklist_bt  = pbt.BedTool(blacklist) 

def exclude_regions_by_bedtool(bt, filter): 
   """ take a bedtool and filter bedtool to exclude any intervals that overlap the blacklist """
   if filter.any_hits(bt): 
      return(False) 
   else: 
      return(True) 

files = [p for p in os.listdir(peakdir) if p.endswith(patt) ]

for f in files: 
   bt = pbt.BedTool(peakdir + f) 
   name = f.split('.')[0]
   bt.merge(d=100)
   bt_filtered = bt.filter(exclude_regions_by_bedtool, filter = blacklist_bt) 
   bt_filtered.saveas('/magnuson-lab/jraab/scratch/tmp.bed')
   df = pd.read_table('/magnuson-lab/jraab/scratch/tmp.bed', names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'x', 'y', 'z'] ) 
   df = df.ix[:,0:6]
   print df.head()
   df['name'] = [name +'_'+ str(i) for i in range(df['name'].shape[0] ) ] 
   df.to_csv( peakdir + 'cleaned/' + name + '_cf.bed', index=False, header=False, sep='\t' ) 



