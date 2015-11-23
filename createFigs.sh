#!/bin/bash 

#Process chip peaks (real peaks already called) 
echo "Preprocessing chip Peaks"
Rscript code/chip_analysis/define_enhancers.R

echo "Annotating Peaks"
Rscript code/chip_analysis/annotate_arid_peaks.R
bash code/chip_analysis/make_bed_of_peaks.sh

echo "Generating Pairwise Correlations for ENCODE DATA"
#Generate pairwise enrichment values for figure creation
Rscript code/chip_analysis/pairwise_enrichments.R

echo " Running Deeptools " 
# DeepTools 
sh code/chip_analysis/arid_coverage_aridpeaks_dt.sh  # Supplemental Figure 4 -1
sh code/chip_analysis/arid_coverage_enhancers_dt.sh  # Figure 3B
sh code/chip_analysis/arid_coverage_genes_dt         # Figure 3A


#Create Figures
# The actual figures in the paper were pasted together in inkscape to make the full panels
# This will dump figures into figures/Fig#_panel.png

" Echo Creating Figures - This may take a while"
Rscript code/figures/fig1_rnaseq_overview.R
Rscript code/figures/fig2_rnaseq_avj.R
Rscript code/figures/fig3_chipseq_ov.R
Rscript code/figures/fig4_chipseq_multi.R
Rscript code/figures/fig5_factor_similarities.R
Rscript code/figures/fig6_metaplots.R
Rscript code/figures/fig7_chip_rna_integrate.R
Rscript code/figures/fig8_rtpcr.R
Rscript code/figures/fig9_id_factors.R

# Extra data tables and things requested by reviewers that did not create a specific figure 
# List of genes that are bound by each arid and the transcription change associated with that arid
# geneName, #_of_peaks_assigned, log2foldchange, p.adjust
Rscript code/rna_chip_analysis/link_chippeaks_toexp_changes.R

# Reviewer asked for amounts of ARIDS in hepg2 relative to other cell lines 
# Used data from encode to estimate this
Rscript code/rna_analysis/calculate_tpms_fromencode.R
Rscript code/rna_chip_combined/all_genes_table_chip_and_rna.R

# Reviewer asked to see chip folloowing knockdown. 
# this is new figure 6
Rscript code/figures/fig6_chipkd.R

# Validate with second set of siRNAs
Rscript code/rtpcr/plot_supp_combo.R
Rscript code/rtpcr/plot_supp_val.R


# Pathway analysis
Rscript code/rna_analysis/pathway_analysis_aloneVJoint.R 
# Output of these files was copied into excel by hand - raw output can be found in output/diffexp/tables/go/

