#!/bin/bash
# use deeptools computematrix and heatmapper to show arid binding
# at arid peaks

bed=output/misc/hepg2_enhancers_fordt.bed

matrix_fdir=output/dt_matfiles/

if [ ! -d $matrix_fdir ]; then
   mkdir -p $matrix_fdir
fi

# do k27ac separatly first to get sorted order
k27sig=data/external/wgEncodeBroadHistoneHepg2H3k27acStdSig.bw
#computeMatrix reference-point \
#    --scoreFileName $k27sig \
#    --regionsFileName $bed \
#    --referencePoint center \
#    --beforeRegionStartLength 2500 \
#    --afterRegionStartLength 2500 \
#    --binSize 10 \
#    --sortRegions no \
#    --numberOfProcessors 4 \
#    --outFileName ${matrix_fdir}enhancers_k27ac_mat
#
#heatmapper \
#    --matrixFile ${matrix_fdir}enhancers_k27ac_mat \
#    --outFileName ${matrix_fdir}enhancers_k27ac.pdf \
#    --sortRegions descend \
#    --colorMap Greys \
#    --whatToShow "heatmap only" \
#    --outFileSortedRegions tmp_k27sorted_enhancers.bed
#

# Next do the arids using the above sorted bed file
signal=$(find data/chip/processed/bw \( -name "arid*.bw" -o -name "snf5*.bw" \) )

#  for i in $signal; do
#     echo $i
#     ofn=$(basename $i |sed 's/.bw//g')
#     computeMatrix reference-point \
#        --scoreFileName $i \
#        --regionsFileName tmp_k27sorted_enhancers.bed \
#        --referencePoint center \
#        --beforeRegionStartLength 2500 \
#        --afterRegionStartLength 2500 \
#        --binSize 10 \
#        --sortRegions no \
#        --numberOfProcessors 4 \
#        --outFileName ${matrix_fdir}enhancers_${ofn}_mat
#  done

mapfiles=$(find $matrix_fdir -name "enhancers_*_mat")
cmap=(Reds Blues Greens Greys Purples)
v=0
echo ${cmap[v]}
for i in $mapfiles; do
   echo $i
   ofn=$(basename $i |sed 's/_mat//g').pdf
   heatmapper \
   --matrixFile $i \
   --outFileName ${matrix_fdir}${ofn} \
   --sortRegions no \
   --colorMap ${cmap[v]} \
   --whatToShow "heatmap only"
   v=$v+1
done

#clean up the tmp file
rm tmp_k27sorted_enhancers.bed
