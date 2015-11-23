#!/bin/bash

cat output/diffexp/tables/expressed_genes_ordered_by_exp.bed > \
    output/diffexp/tables/combined_genes.bed
cat output/diffexp/tables/not_expressed_genes.bed >> \
    output/diffexp/tables/combined_genes.bed

bed=output/diffexp/tables/combined_genes.bed
pdir=data/chip/processed/bw/

bwfiles=$(find $pdir \( -name "arid*.bw" -o -name "snf5*.bw" \) )

echo $bwfiles
odir=output/dt_matfiles/ # change this to be same as on osx
if [ ! -d $odir ] ; then
    mkdir =p $odir
fi

for i in $bwfiles; do
    echo $i
    ofn=$(basename $i | sed  's/.bw//g')
    computeMatrix scale-regions \
        --scoreFileName $i \
        --regionsFileName $bed \
        --regionBodyLength 3000 \
        --startLabel TSS \
        --endLabel TES \
        --beforeRegionStartLength 2000 \
        --afterRegionStartLength 1000 \
        --numberOfProcessors 3 \
        --skipZeros \
        --outFileName ${odir}genes_${ofn}_mat
done


mapfiles=$(find ${odir} -name "genes_*_mat")
echo $mapfiles
cmap=(Reds Blues Greens Purples)
v=0

echo ${cmap[0]}

for i in $mapfiles; do
    echo $i
    heatmapper \
    --matrixFile $i \
    --outFileName ${i}.pdf \
    --sortRegions  no \
    --colorMap ${cmap[v]} \
    --whatToShow "heatmap only"
    v=$v+1

done
