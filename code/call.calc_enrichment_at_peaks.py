#!/usr/bin/env python
import subprocess
import os
import re
import time
import yaml

# wrapper script to repeatedly call my enrichment mapper
# The directories are not set up to be portable - but the only script this depends
# on is submitCoverageFrame.sh which calls coverageFrame.py
# These are provided in the utils/ directory

peakdir = '/magnuson-lab/jraab/analysis/swi_snf_final/output/macs_peaks/cleaned/'

def allEncodeOverPeaks(peaks, peakdir):
   bamdir = '/magnuson-lab/jraab/ENCODE/datafiles/HepG2/combined/'
   bams = [b for b in os.listdir(bamdir)if b.endswith('.bam') ]
   out = '/magnuson-lab/jraab/analysis/swi_snf_final/output/encode_coverages/'

   if not os.path.exists(out):
      os.mkdir(out)
   yf = open('/magnuson-lab/jraab/ENCODE/datafiles/HepG2/hepg2_input_index.yaml')
   yobj = yaml.load(yf)
   for input, p in yobj.iteritems():
      for ip in p:
         ip_f = bamdir+ip
         print input, ip
         input_f = bamdir+input
         for s in peaks:
            bamname = ip.split('.')[0].split('Hepg2')[1]
            pname = s.split('_')[0]
            print bamname, pname
            qsub_cmd = 'qsub -V -v P='+peakdir+s+',B='+ip_f+',I='+input_f+',BNAME='+bamname+',O='+out+',PNAME='+pname+' submitCoverageFrame.sh'
            print qsub_cmd.split()
            subprocess.call(qsub_cmd.split(), shell=False)
            time.sleep(3)
         time.sleep(10)


   yf = open('/magnuson-lab/jraab/ENCODE/datafiles/HepG2/haib.yml')
   yobj = yaml.load(yf)
   for input, p in yobj.iteritems():
      for ip in p:
         ip_f = bamdir+ip
         print input, ip
         input_f = bamdir+input
         for s in peaks:
            bamname = ip.split('.')[0].split('Hepg2')[1]
            pname = s.split('_')[0]
            print bamname, pname
            qsub_cmd = 'qsub -V -v P='+peakdir+s+',B='+ip_f+',I='+input_f+',BNAME='+bamname+',O='+out+',PNAME='+pname+' submitCoverageFrame.sh'
            print qsub_cmd.split()
            subprocess.call(qsub_cmd.split(), shell=False)
            time.sleep(3)
         time.sleep(10)

#pass in a list of peak summits and a directory they came ffrom

peakdir = '/magnuson-lab/jraab/analysis/swi_snf_final/output/macs_peaks/cleaned/'
summits = [ s for s in os.listdir(peakdir) if s.endswith('_cf.bed')]
print summits
#allEncodeOverPeaks(summits, peakdir)

multi_summit = [s for s in os.listdir(peakdir) if s.endswith('multi_bound_peaks.csv') ]
#allEncodeOverPeaks(multi_summit, peakdir)
# do the dnase separately
summits.append(multi_summit[0])
for peaks in summits:
   bamdir = '/magnuson-lab/jraab/ENCODE/datafiles/HepG2/combined/'
   out = '/magnuson-lab/jraab/analysis/swi_snf_final/output/encode_coverages/'
   signal_file = 'wgEncodeOpenChromDnaseHepg2Aln.sorted.merged.bam'
   bamname = 'DNase'
   pname = peaks.split('_')[0]
   ip_f = bamdir+signal_file
   qsub_cmd = 'qsub -V -v P='+peakdir+peaks+',B='+ip_f+',I=NULL,BNAME='+bamname+',OUT='+out+',PNAME='+pname+' submitCoverageFrame.sh'
   print qsub_cmd.split()
   subprocess.call(qsub_cmd.split(), shell = False)
   time.sleep(3)
