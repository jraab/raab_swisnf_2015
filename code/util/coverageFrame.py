#!/usr/bin/env python
#########################
# Given a bed file of interval info and a bamfile of reads
# calculate the enrichment over each interval - if named, intervals will
# retain name, else get a random name
# I usually use a wrapper script to call this using a cluster
# see call.calc_enrichment_at_peaks.py
# This was used to generate files in
#    output/encode_coverages/
#####################################

import HTSeq
import pandas as pd
import numpy as np
import argparse
import sys
import os
import gzip
import pysam
import tables
import re
import hashlib
##############################
#arguments
parser= argparse.ArgumentParser(description='calculate enrichment at a set of intervals')
parser.add_argument('-p', help='peakfile')
parser.add_argument('-b', help='bamfile')
parser.add_argument('-i', nargs='?', help='input file to normalize enrichment to')
parser.add_argument('-pname', help='peakfile name - for naming output')
parser.add_argument('-bamname', help='bamfile name - for naming output')
parser.add_argument('-pkwin', default=250, help='size of enrichment window')
parser.add_argument('-o', nargs='?', default='.',  help='output directory')
parser.add_argument('-w', default = 5000, help='bp to output in coverage file - will be binned by -b')
parser.add_argument('-bg', default = 20000, help='size of the background window to calculate enrichment')
parser.add_argument('-m', action='store_true', help='calculate midpoint of interval? to center signal')
parser.add_argument('-binsize', default = 50, help='size of bins for output of coverage')
parser.add_argument('-binaction', choices=['mean', 'median', 'sum'], default = 'mean', help='median, mean, or sum for calculating read coverage within a bin')
args = parser.parse_args()
#############################################
args.pkwin = int(args.pkwin)
args.w = int(args.w)
args.binsize = int(args.binsize)
args.bg = int(args.bg)
if not os.path.exists(args.o):
   os.makedirs(args.o)

#function defs

def getsep(fname):
   """
   Check if file is tab or comma separated - doesn't currently support other delimiters'
   """
   with open(fname) as fn:
      lines = 0
      if lines < 10:
         for line in fn:
            if line.startswith('#'):
               continue
            x = line.strip().split('\t')
            if len(x) == 1:
               try:
                  x = line.strip().split(',')
                  if len(x) > 1:
                     return(',')
                  else:
                     sys.exit('No valid field separator shown')
               except:
                  sys.exit('No valid field separator shown')
            elif len(x) > 1:
               return('\t')
         lines += 1

def readGenes (fname):
   """
   Take a file name that is bed-like (needs at least 3 columns) - will conserve strand and name if they exist
   Name will be the dictionary key for the output coverages
   """
   genes = {}
   i = 1
   fieldsep = getsep(fname)
   patt = re.compile(r'^chr$')
   with open(fname) as fn:
      for line in fn:
         if line.startswith('#'):
            continue
         try:
            chrom,start,end,name,score,strand = line.strip().split(fieldsep)[:6]
         except:
            try:
               a = line.strip().split(fieldsep)
               length = len(a)
               ext = 6-length
               a.extend(['.'] * ext)
               chrom, start, end, name, score, strand = a
            except:
               sys.error(fname, 'Not a valid bed')
         if name == '.':
            name = 'peak'+str(i)
         else:
            name = os.path.basename(name) #macs puts the full path in this name field, if it hasn't been cleaned up yet this will
         if strand not in ('+', '-', '.'):
            strand = '.'
         if args.m:
            try:
               start = (int(start)+int(end)) / 2
               end = int(start +1 )
            except:
               continue
         else:
            #this should be used to deal with TSS
            start = int(start)
            end = int(start)+1
         iv = HTSeq.GenomicInterval(chrom.replace('"',''), int(start), int(end), strand)
         genes[str(name)] = iv
         i += 1
   return(genes)

def bamCoverage (pos_dic, bamfile, halfwinwidth, fragmentsize = 200 ):
   frame = {}
   bf = pysam.Samfile(bamfile, 'rb', check_header=False, check_sq=False)
   for key in pos_dic:
      p = pos_dic[key]
      profile = np.zeros(2*halfwinwidth, dtype='float')
      window = HTSeq.GenomicInterval(p.chrom, p.start_d-halfwinwidth-fragmentsize, p.start_d+halfwinwidth+fragmentsize, '.')
      if window.start < 0:
         print 'window issue'
         continue
      for a in bf.fetch(p.chrom, window.start, window.end):
         if a.is_reverse:
            strand = '-'
         else:
            strand = '+'
         almnt = HTSeq.GenomicInterval(window.chrom, a.pos, a.pos+fragmentsize, strand)
         if p.strand == '+' or p.strand == '.':
            start_in_window = almnt.start - p.start_d + halfwinwidth
            end_in_window = almnt.end - p.start_d + halfwinwidth
         else:
            start_in_window = p.start_d + halfwinwidth - almnt.end
            end_in_window = p.start_d + halfwinwidth - almnt.start
         start_in_window = max(start_in_window, 0)
         end_in_window = min(end_in_window, 2*halfwinwidth)
         if start_in_window >= 2*halfwinwidth or end_in_window < 0 :
            continue
         profile[start_in_window:end_in_window] += 1
      frame[key] = profile
   return(frame)



def excludeRegions(genedict):
   """
   take a dictionary of gene-intervals and check them against known blacklisted regions - returning only those
   that do not overlap the restricted areas
   """
   e1 = gzip.open('/magnuson-lab/jraab/annotations/wgEncodeDukeMapabilityRegionsExcludable.bed.gz')
   e2 = gzip.open('/magnuson-lab/jraab/annotations/wgEncodeDacMapabilityConsensusExcludable.bed.gz')

   excludedRegions = []
   cleanGenes ={}
   for l1 in e1:
      chr,start,end,name,score,strand=l1.split('\t')
      interval = HTSeq.GenomicInterval(chr, int(start), int(end), '.')
      excludedRegions.append(interval)
   for l2 in e2:
      chr,start,end,name,score,strand=l2.split()
      interval = HTSeq.GenomicInterval(chr, int(start), int(end), '.')
      excludedRegions.append(interval)
   for key in genedict:
      gene = genedict[key]
      hit = 0
      for interval in excludedRegions:
         if interval.overlaps(gene):
            hit=1
      if hit==0:
         cleanGenes[key] = gene
   return(cleanGenes)


def readDepth(bamfile):
   try:
      bf = pysam.Samfile(bamfile, 'rb', check_header=False, check_sq=False)
   except:
      print bamfile
      sys.exit()
   return(bf.mapped)

def rpm(val, reads):
   return( val * 1e6/reads)

def enrichmentScore(numpy_vec, peak_window):
   """
   Given a numpy array calculate how much enrichment over the background window is
   found in the peak_window
   Simplistic socring of log2(median(small_window)+1/median(background)+1)
   """
   mid = int(len(numpy_vec)/2)
   subset = numpy_vec[mid-peak_window: mid+peak_window]
   num = float((np.mean(subset)))
   den = float((np.mean(numpy_vec)))
   enrichment = np.log2(num/den)
   return(enrichment)

def subset_vec(numpy_vec, winsize):
   midpt = int(len(numpy_vec)/2)
   halfwin = int(winsize/2)
   try:
      return(numpy_vec[midpt-halfwin:midpt+halfwin])
   except:
      print midpt, type(midpt), halfwin, type(halfwin)

## this section is not used
#def saveHDF(from_dict, fname, group):
#   length_vals = len(from_dict.values()[0])
#   class Coverage(tables.IsDescription):
#      name = tables.StringCol(24)
#      vals = tables.Float64Col(shape=(1,length_vals) )
#   f = tables.open_file(fname, 'w')
#   t = f.create_table(f.root, 'coverages', Coverage)
#   cov = t.row
#   for k,v in from_dict.iteritems():
#      cov['name'] = k
#      cov['vals'] = v
#      cov.append()
#   t.flush()
#   f.root.coverages.cols.name.create_index()
#   f.close()

##############################################################
def main():
   genes = readGenes(args.p)
   genes = excludeRegions(genes)
   peak_depth  = readDepth(args.b)
   peak_coverage = bamCoverage(genes, args.b, int(args.bg)/2)
   peak_coverage = {k:v+1. for k,v in peak_coverage.iteritems()}
   peak_coverage = {k:rpm(v, peak_depth) for k,v in peak_coverage.iteritems() }
   # get info for the input file if included
   # There is some ugly error handling going on below related to checkedin
   # for keys in the input and ip dictionaries these should be identical
   #
   if args.i:
      input_depth = readDepth(args.i)
      input_coverage = bamCoverage(genes, args.i, int(args.bg)/2)
      input_coverage = {k:v+1. for k,v in input_coverage.iteritems() }
      try:
         input_coverage = {k:rpm(v, input_depth) for k,v in input_coverage.iteritems() }
      except:
         'input_coverage failed'
   #normalize to input (as IP+1/Input+1)

   output_coverage = {}
   for k,v in peak_coverage.iteritems():
      if args.i and k in input_coverage.keys():
         try:
            input = input_coverage[k]
            out_val = v / input  # this was failing in previous version becaues i had called float(v)/float(input)
            # seems it was unnecessary
         except:
            out_val = v
            print 'not normalized'
         output_coverage[k] = out_val
      elif args.i:
         print 'key not in input, dropping %s ' % k
      else:
         output_coverage[k] = v
   peak_coverage = None
   input_coverage = None
   #calculate enrichment score over TSS
   enrichments = {k:enrichmentScore(v, int(args.pkwin)) for k,v in output_coverage.iteritems() }
   subsets = {k:subset_vec(v,args.w) for k,v in output_coverage.iteritems()}
   #convert to dataframe for saving
   enrich_df = pd.DataFrame.from_dict(enrichments, orient='index')
   enrich_df.index = [x.replace('"', '') for x in enrich_df.index.values ]
   enrich_df.index.name = 'gene'
   enrich_df.columns = ['enrichment']
   cov_dir = args.o+'/coverages/'
   enrich_dir = args.o+'/enrichments/'

   for a in cov_dir,enrich_dir:
      if not os.path.exists(a):
         os.mkdir(a)

   out_fn = cov_dir+args.pname+'_p_'+ args.bamname+'_s_coverage.h5'
   enrich_fn = enrich_dir+args.pname+'_p_'+args.bamname+'_s_enrichment.csv'
   subsets = pd.DataFrame.from_dict(subsets, orient='index')
   subsets.index = [x.replace('"', '') for x in subsets.index.values ]
   subsets.index.name = 'gene'
   #trim the data set by the bin amount to cut down on space on disk these take up
   if args.binaction == 'mean':
      print 'mean'
      binned = subsets.transpose().groupby(lambda x: x/args.binsize).mean().transpose()
   elif args.binaction == 'median':
      print 'median'
      binned = subsets.transpose().groupby(lambda x: x/args.binsize).median().transpose()
   elif args.binaction == 'sum':
      print 'sum'
      binned = subsets.transpose().groupby(lambda x: x/args.binsize).sum().tranpsose()
   binned.to_hdf(out_fn, 'coverages', mode='w' )
   #store = pd.HDFStore(out_fn)
   #store['coverages'] = binned
   #store.close()
   enrich_df.to_csv(enrich_fn, sep=',')

if __name__ == '__main__':
   main()
