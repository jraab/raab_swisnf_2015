#!/bin/bash
# use deeptools computematrix and heatmapper to show arid binding
# at arid peaks

bed=output/macs_peaks/cleaned/all_peaks_fordt.bed
signal=$(find data/chip/processed/bw \( -name "arid*.bw" -o -name "snf5*.bw" \) )
matrix_fdir=output/dt_matfiles/
if [ ! -d $matrix_fdir ]; then
   mkdir -p $matrix_fdir
fi

for i in $signal; do
   echo $i
   ofn=$(basename $i |sed 's/.bw//g')
   computeMatrix reference-point \
      --scoreFileName $i \
      --regionsFileName $bed \
      --referencePoint center \
      --beforeRegionStartLength 2500 \
      --afterRegionStartLength 2500 \
      --binSize 10 \
      --sortRegions no \
      --numberOfProcessors 2 \
      --outFileName ${matrix_fdir}aridpeaks_${ofn}_mat
done

mapfiles=$(find $matrix_fdir -name "aridpeaks_*_mat")
cmap=(Reds Blues Greens Purples)
v=0
for i in $mapfiles; do
   echo $i
   ofn=$(basename $i |sed 's/_mat//g').pdf
   heatmapper \
   --matrixFile $i \
   --outFileName ${matrix_fdir}${ofn} \
   --sortRegions no \
   --colorMap ${cmap[v]} \
   --whatToShow "plot and heatmap"
   v=$v+1
done
