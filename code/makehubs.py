#!/usr/bin/env python
import os
from trackhub import Track, default_hub, CompositeTrack, ViewTrack
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
uploadbase = "http://trackhubs.its.unc.edu/magnuslb/jraab/swi_snf/bw/"
remote_dir = '/proj/magnuslb/genome_browser/jraab/swi_snf/bw/'
datadir = '/magnuson-lab/jraab/swi_snf_final/data/'

#setup dict to map antibody to colors
maps = {'input': '105,105,105','arid1a':'228,26,28', 'arid1b': '55,126,184', 'arid2':'77,175,74', 'snf5': '152,78,163' } 


def get_name_and_rep(string): 
   rep = string.split('_')[1].split('.')[0]
   name = string.split('_')[0]
   return(name,rep)

#set up composite track to hold chipseq info 
comp = CompositeTrack(
   name='swi_snf_chipseq',
   short_label='swi_snf_super',
   long_label='swi_snf_chipseq supertrack',
   tracktype='bigWig',
   dragAndDrop='subtracks',
   visibility='full',
   viewLimits='1:10',
   autoScale='off',
   maxHeightPixels='100:30:10',
   priority=1
   )


#######################
#make the tracks and add to the appropriate groups during loo
chipseq_bwfiles = [f for f in os.listdir(datadir+'/chip/processed/bw/') if f.endswith('.bw')] 
chipseq_bbfiles = [b for b in os.listdir(datadir+'/chip/processed/bb/') if b.endswith('.bb')]
signal_view = ViewTrack(
      name = 'Chipseq_signal', 
      view = 'Signal', 
      visibility = 'full', 
      tracktype = 'bigWig', 
      short_label = 'signal')
bed_view = ViewTrack(
      name = 'Chipseq_bed', 
      view = 'Bed', 
      visibility = 'dense',
      tracktype = 'bigBed',
      short_label = 'peaks')
comp.add_view(signal_view)
comp.add_view(bed_view)


for f in chipseq_bwfiles: 
   name, rep = get_name_and_rep(f)
   bname = os.path.basename(f) 
   col = maps[name]  
   track = Track(
      name = '%s_%s_signal' % (name,rep) , 
      tracktype='bigWig',
      local_fn = datadir + f, 
      remote_fn = remote_dir + bname, 
      url = uploadbase+'bw/'+ bname,
      shortLabel = '%s_%s' %( name,rep) , 
      color = col)
   signal_view.add_tracks(track)

for b in chipseq_bbfiles:  
   bname = os.path.basename(b)
   name = bname.split('.')[0]
   col = maps[name]
   track = Track(
      name = '%s_peak' % name,
      tracktype='bigBed',
      local_fn = datadir+ f, 
      remote_fn = remote_dir + bname, 
      url = uploadbase + 'bb/'+bname,  
      shortLabel = '%s_peaks' % name, 
      color=col) 
   bed_view.add_tracks(track)

#need to add in the RNAseq to same composite groups
trackdb.add_tracks(comp)
print trackdb
hub.render()

kwargs = dict(host='kure.its.unc.edu', user='jraab') 
upload_hub(hub=hub, **kwargs) 
for track, level in hub.leaves(Track):  
   upload_track(track=track, **kwargs)
