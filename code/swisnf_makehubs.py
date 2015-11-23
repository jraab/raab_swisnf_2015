#!/usr/bin/env python
import os
from trackhub import Track, default_hub, CompositeTrack, ViewTrack, SuperTrack
import re
from trackhub.upload import upload_hub, upload_track

# important directories
hub, genomes_file, genome, trackdb = default_hub(
   hub_name = 'swi_snf', 
   genome = 'hg19', 
   short_label = 'swi_snf', 
   long_label = 'swi_snf_genomics', 
   email = 'jesse.r.raab@gmail.com' )
hub.url = 'http://trackhubs.its.unc.edu/magnuslb/jraab/swi_snf/hub.txt'
hub.remote_fn = '/proj/magnuslb/genome_browser/jraab/swi_snf/hub.txt'
uploadbase = "http://trackhubs.its.unc.edu/magnuslb/jraab/swi_snf/"
remote_dir = '/proj/magnuslb/genome_browser/jraab/swi_snf/bw/'
datadir = '/magnuson-lab/jraab/analysis/swi_snf_final/data/'

#setup dict to map antibody to colors chipmaps = {'input': '105,105,105',
chipmaps = { 'arid1a':'228,26,28', 
             'arid1b': '55,126,184', 
             'arid2':'77,175,74', 
<<<<<<< HEAD
             'snf5': '152,78,163', 
             'input': '99,99,99'} 

rnamaps = { 'NTG_plus': '99,99,99', 
            'NTG_minus':'189,189,189',
            'ARID1A_plus': '222,45,38' , 
            'ARID1A_minus':'252,146,114',
            'ARID1B_plus': '49,130,189', 
            'ARID1B_minus': '158,202,225', 
            'ARID2_plus': '49,163,84',
            'ARID2_minus': '161,217,155' } 
# Top level group - Supertrack 
super = SuperTrack(
   name = 'swi_snf', 
   long_label = 'Swi Snf Supertrack', 
   short_label = 'snf_snf_supertrack'
   )

#set up composite track to hold chipseq info 
comp_chip = CompositeTrack(
   name='swi_snf_chipseq',
   short_label='swi_snf_chipseq',
   long_label='swi_snf_chipseq', 
   tracktype = 'bigWig',
   dragAndDrop='subtracks', 
   viewLimits='4:30',
   autoScale='off',
   maxHeightPixels='100:15:10',
   priority=2   
   )

comp_rna = CompositeTrack(
   name = 'swi_snf_rnaseq', 
   short_label = 'swi_snf_rnaseq', 
   long_label = 'swi_snf_rnaseq siRNA', 
   tracktype = 'bigWig',
   dragAndDrop = 'subtracks', 
   visibility = 'full', 
   viewLimits = '0:5', 
   autoScale = 'off', 
   maxHeightPixels = '100:20:10', 
   priority=1 
   )  

# Add the compositive views ot the super track 
# ----------------------------------------------------------------------------
super.add_track(comp_chip)
super.add_track(comp_rna) 
=======
             'snf5': '152,78,163' } 

#rnamaps = { 'NTG_plus': '99,99,99', 
#            'NTG_minus':'189,189,189',
#            'ARID1A_plus': '222,45,38' , 
#            'ARID1A_minus':'252,146,114',
#            'ARID1B_plus': '49,130,189', 
#            'ARID1B_minus': '158,202,225', 
#            'ARID2_plus': '49,163,84',
#            'ARID2_minus': '161,217,155' } 
#
#set up composite track to hold chipseq info 
comp = CompositeTrack(
   name='swi_snf_chipseq',
   short_label='swi_snf',
   long_label='swi_snf_chipseq', 
   tracktype = 'bigWig',
   dragAndDrop='subtracks',
   visibility='full',
   viewLimits='3:40',
   autoScale='off',
   maxHeightPixels='100:20:10',
   priority=2   
   )

#comp_rna = CompositeTrack(
#   name = 'swi_snf_rnaseq', 
#   short_label = 'swi_snf_rnaseq', 
#   long_label = 'swi_snf_rnaseq siRNA', 
#   dragAndDrop = 'subtracks', 
#   visibility = 'full', 
#   viewLimits = '0:5', 
#   autoScale = 'off', 
#   maxHeightPixels = '100:20:10', 
#   priority = 1, 
#   
#   )  

# Add the compositive views ot the super track 
# ----------------------------------------------------------------------------

>>>>>>> e33a5d89c0f5d75afae01dde9452063e2287f8d6


# Creat the views and add views to composite appropriate composite tracks
# -----------------------------------------------------------------------------
#make the tracks and add to the appropriate groups during loo
chipseq_bwfiles = [f for f in os.listdir(datadir+'/chip/processed/bw/') if f.endswith('.bw')] 
chipseq_bbfiles = [b for b in os.listdir(datadir+'/chip/processed/bb/') if b.endswith('.bb')]
<<<<<<< HEAD
rnaseq_bwfiles  = [r for r in os.listdir(datadir+'/rna/processed/bw/') if r.endswith('.bw')] 

print chipseq_bwfiles
print chipseq_bbfiles
print rnaseq_bwfiles
signal_view = ViewTrack(
=======
#rnaseq_bwfiles  = [r for r in os.listdir(datadir+'/rna/processed/bw/') if r.endswith('.bw')] 
chip_signal_view = ViewTrack(
>>>>>>> e33a5d89c0f5d75afae01dde9452063e2287f8d6
      name = 'Chipseq_signal', 
      view = 'Signal', 
      visibility = 'full', 
      tracktype = 'bigWig', 
      short_label = 'signal')
<<<<<<< HEAD
bed_view = ViewTrack(
=======
chip_bed_view = ViewTrack(
>>>>>>> e33a5d89c0f5d75afae01dde9452063e2287f8d6
      name = 'Chipseq_bed', 
      view = 'Bed', 
      visibility = 'dense',
      tracktype = 'bigBed',
      short_label = 'peaks')
<<<<<<< HEAD

rna_view = ViewTrack(
      name = 'RNAseq_signal', 
      view = 'RNA', 
      visibility = 'full', 
      tracktype = 'bigWig', 
      short_label = 'Rnaseq_signal') 

comp_chip.add_view(signal_view)
comp_chip.add_view(bed_view)
comp_rna.add_view(rna_view)
=======
#rna_signal_view = ViewTrack(
#      name = 'RNAseq_signal', 
#      view = 'Signal', 
#      visibility = 'full', 
#      tracktype = 'bigwig', 
#      short_label = 'signal') 
#
comp.add_view(chip_signal_view)
comp.add_view(chip_bed_view)
#comp_rna.add_view(rna_signal_view)
>>>>>>> e33a5d89c0f5d75afae01dde9452063e2287f8d6

# Loop through all files and create tracks - adding to appropriate view
# -------------------------------------------------------------------------------

for f in chipseq_bwfiles:  
   bname = os.path.basename(f) 
   name = bname.split('.')[0]
   col = chipmaps[name]  
   track = Track(
      name = '%s_signal' % (name) , 
      tracktype='bigWig',
      local_fn = datadir + 'chip/processed/bw/'+ f, 
      remote_fn = remote_dir + bname, 
      url = uploadbase + 'bw/' + bname,
      shortLabel = '%s' %(name) , 
      color = col) 
<<<<<<< HEAD
   signal_view.add_tracks(track)
=======
   chip_signal_view.add_tracks(track)
>>>>>>> e33a5d89c0f5d75afae01dde9452063e2287f8d6

for b in chipseq_bbfiles:  
   bname = os.path.basename(b)
   name = bname.split('.')[0]
   col = chipmaps[name]
   track = Track(
      name = '%s_peak' % name,
      tracktype='bigBed',
      local_fn = datadir + 'chip/processed/bb/' + f, 
      remote_fn = remote_dir + bname, 
      url = uploadbase + 'bb/'+bname,  
      shortLabel = '%s_peaks' % name, 
      color = col) 
<<<<<<< HEAD
   bed_view.add_tracks(track)

for r in rnaseq_bwfiles: 
   bname = os.path.basename(r)
   name, strand = bname.split('.')[0].split('_')
   col = rnamaps[bname.split('.')[0]] 
   track = Track(
      name = '%s_%s_signal' % (name, strand), 
      tracktype = 'bigWig', 
      local_fn = datadir + 'rna/processed/bw/' + r, 
      remote_fn = remote_dir + bname, 
      url = uploadbase + 'bw/' + bname,
      shortLabel = '%s_%s_rnasignal' % (name, strand) , 
      color = col) 
   rna_view.add_tracks(track) 

# -----------------------------------------------------------------------------
# Add everyhing in the supertrack to the db and upload data to Kure
trackdb.add_tracks(super)

print trackdb
hub.render()

kwargs = dict(host='kure.its.unc.edu', user='jraab')
=======
   chip_bed_view.add_tracks(track)

#for r in rnaseq_bwfiles: 
#   bname = os.path.basename(r)
#   name, strand = bname.split('.')[0].split('_')
#   col = rnamaps[bname.split('.')[0]] 
#   track = Track(
#      name = '%s_%s_signal' % (name, strand), 
#      tracktype = 'bigWig', 
#      local_fn = datadir + 'rna/processed/bw/' + bname, 
#      remote_fn = remote_dir + bname, 
#      url = uploadbase + 'bw/' + bname,
#      shortLabel = '%s_%s_rnasignal' % (name, strand) , 
#      color = col) 
#   rna_signal_view.add_tracks(track) 
#
# -----------------------------------------------------------------------------
# Add everyhing in the supertrack to the db and upload data to Kure

trackdb.add_tracks(comp)
print trackdb
hub.render()

kwargs = dict(host='kure.its.unc.edu', user='jraab') 
>>>>>>> e33a5d89c0f5d75afae01dde9452063e2287f8d6
upload_hub(hub=hub, **kwargs) 
for track, level in hub.leaves(Track):  
   upload_track(track=track, **kwargs)
