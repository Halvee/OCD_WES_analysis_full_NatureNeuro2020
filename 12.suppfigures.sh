#!/bin/bash

# needed for some python scripts - access to installed packages
export PYTHONPATH="$PWD/lib/python/:$PYTHONPATH"

# needed for some R scripts - access to installed packages
export R_LIBS=$PWD/lib/R/

# make sure output dir exists
mkdir -p results/suppfigures/

# figure S1 : UMAP of input case/control PCA for collapsing
convert -density 300 -quality 100 \
data/suppfigures/FigureS1.sample_UMAP.pdf \
results/suppfigures/FigureS1.sample_UMAP.jpg

# figure S2 : images in IGV of alns underlying case SLITRK5 LoF/MisD calls
for num in "1" "2" "3" "4" "5" "6"
do
  convert -density 300 -quality 100 \
  data/suppfigures/SLITRK5_IGV_plots_20210420_pdfs/SLITRK5_IGV_plots_20210420.p${num}.pdf \
  results/suppfigures/FigureS2.SLITRK5_IGV.${num}.png

done

# figure S3 : 'lollipop' plot of SLITRK5 results
convert -density 300 -quality 100 \
data/suppfigures/SLITRK5_lollipop_plot.formatted.pdf \
results/suppfigures/FigureS3.SLITRK5_lollipop_plot.png

# figure S4 : PCA of European ancestry samples used for LoF rate comparisons :
#             all cases from cluster 0, cluster 3 or cluster 4;
#             the subset of controls from these same 
#             clusters that have phenotype listing of 
#             'control' or 'healthy family member'
bin/Rscript src/suppfigures/collapsing_pca_plots_single.R \
data/clps_rate/OCD_clps_rate.allcases.evec \
PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8 \
results/suppfigures/FigureS4.clps_caco_rate_pca_PC1-PC8.pdf
convert -density 300 -quality 100 \
results/suppfigures/FigureS4.clps_caco_rate_pca_PC1-PC8.pdf \
results/suppfigures/FigureS4.clps_caco_rate_pca_PC1-PC8.jpg

# figure S5 : case/control rate tests 
# (LOEUF, linear regression, LoF/synon fisher's exact test)
rsync \
results/clps_rate/OCD_clps_rate.allcases.LOF_LOEUF.linear_fet.pdf \
results/suppfigures/FigureS5.LOF_LOEUF_linear_fet.pdf
convert -density 300 -quality 100 \
results/suppfigures/FigureS5.LOF_LOEUF_linear_fet.pdf \
results/suppfigures/FigureS5.LOF_LOEUF_linear_fet.png

# figure S6 : case/control rate tests 
# (pLI, pLI+LOEUF, linear+logistic reg, LoF/synon fisher's exact test)
rsync \
results/clps_rate/OCD_clps_rate.allcases.pLI.pdf \
results/suppfigures/FigureS6.LOF_PLI_LogisticLinearFet.pdf
convert -density 300 -quality 100 \
results/suppfigures/FigureS6.LOF_PLI_LogisticLinearFet.pdf \
results/suppfigures/FigureS6.LOF_PLI_LogisticLinearFet.png

# figure S7 : dnm case/control rate differences across the exome (top),
#             dnm case/control rate diffs of LOF in LOEUF deciles (bottom)
rsync \
results/dnm_rate_caco/OCDfams_2019.cohort.trios.10x.caco_rate.pdf \
results/suppfigures/FigureS7.DNM_caco_rate.pdf
convert -density 300 -quality 100 \
results/suppfigures/FigureS7.DNM_caco_rate.pdf \
results/suppfigures/FigureS7.DNM_caco_rate.png

# figure S8 : dnm obs/exp rate differences across the exome (top),
#             dnm obs/exp rate of  LOF in LOEUF deciles (bottom)
#             expectation based on sequence context
rsync \
results/dnm_rate_caco/OCDfams_2019.cohort.trios.10x.rate.pdf \
results/suppfigures/FigureS8.DNM_rate.pdf
convert -density 300 -quality 100 \
results/suppfigures/FigureS8.DNM_rate.pdf \
results/suppfigures/FigureS8.DNM_rate.png

# figure S9 : dnm rate analysis for LOF variation in LOEUF decile 1 and 2-10,
# partitioned by variant type
# (top: nonsense/splice SNVs, bottom: frameshift indels)
rsync \
results/dnm_rate_caco/OCDfams_2019.cohort.trios.10x.loeuf_bins.non_splice_fs.rate.pdf \
results/suppfigures/FigureS9.DNM_rate.non_splice_frameshift.LOFintoltol.pdf
convert -density 300 -quality 100 \
results/suppfigures/FigureS9.DNM_rate.non_splice_frameshift.LOFintoltol.pdf \
results/suppfigures/FigureS9.DNM_rate.non_splice_frameshift.LOFintoltol.png

# figure S10 : observed/expected de novo mutation rates within
# neurodevelopmental and TS genesets
rsync \
results/dnm_rate_caco/OCDfams_2019.cohort.trios.10x.genesets.misdlof.rate.pdf \
results/suppfigures/FigureS10.DNM_misDLOF_rate_genesets.pdf
convert -density 300 -quality 100 \
results/suppfigures/FigureS10.DNM_misDLOF_rate_genesets.pdf \
results/suppfigures/FigureS10.DNM_misDLOF_rate_genesets.png

# figure S11 : case/control burden of misD DNMs, across LOEUF bins
rsync \
results/dnm_rate_caco/OCDfams_2019.cohort.trios.10x.loeuf_bins.misd.caco_caco_rate.pdf \
results/suppfigures/FigureS11.DNM_misD_LOEUF_bins.pdf 
convert -density 300 -quality 100 \
results/suppfigures/FigureS11.DNM_misD_LOEUF_bins.pdf \
results/suppfigures/FigureS11.DNM_misD_LOEUF_bins.png

# figure S12 : case/control burden of private misD DNMs, seperated
# by MPC bin (0-1, 1-2, >2)
rsync \
results/dnm_rate_caco/OCDfams_2019.cohort.trios.10x.exome.misD_MPC_bins.caco.pdf \
results/suppfigures/FigureS12.misD_MPC_bins_caco.pdf
convert -density 300 -quality 100 \
results/suppfigures/FigureS12.misD_MPC_bins_caco.pdf \
results/suppfigures/FigureS12.misD_MPC_bins_caco.png

# figure S13 : results of tests for assocation between carrier status
# of damaging DNM (LOF/LOEUF<10% OR misD/MPC>2) and key binary features
# (male sex, tics, skin-picking, trichotillomania)
rsync \
results/clinical/OCD.clinical.features.cmp.pdf \
results/suppfigures/FigureS13.DNM_LOFloeuf1_misDmpcgt2_carriers_features_assoc.pdf
convert -density 300 -quality 100 \
results/suppfigures/FigureS13.DNM_LOFloeuf1_misDmpcgt2_carriers_features_assoc.pdf \
results/suppfigures/FigureS13.DNM_LOFloeuf1_misDmpcgt2_carriers_features_assoc.png

# figure S14 : results of LoF rate across LOEUF deciles in
# non-trio male versus non-trio female OCD cases (odds ratios)
rsync \
results/clinical/ocd.nontriocases.malevsfemale.lofrate.LOF_LOEUF.logistic.pdf \
results/suppfigures/FigureS14.LoF_LOEUFdeciles_ORs.malevsfemale.nontrio.pdf
convert -density 300 -quality 100 \
results/suppfigures/FigureS14.LoF_LOEUFdeciles_ORs.malevsfemale.nontrio.pdf \
results/suppfigures/FigureS14.LoF_LOEUFdeciles_ORs.malevsfemale.nontrio.png

# figure S15 : igv, sanger validation for CHD8
convert -density 300 -quality 100 \
data/suppfigures/CHD8_LoF_IGV_Sanger.png \
results/suppfigures/FigureS15.CHD8_LoF_IGV_Sanger.png

# figure S16 : igv, sanger validation for HDAC4
convert -density 300 -quality 100 \
data/suppfigures/HDAC4_LoF_IGV_Sanger.png \
results/suppfigures/FigureS16.HDAC4_LoF_IGV_Sanger.png

# combine all figures into supplemental figures docx
bin/python src/suppfigures/make_supp_figures_docx.py \
--header "Table of Contents" \
--captions-txt data/suppfigures/supp_figure_captions.txt \
--out-docx results/suppfigures/Supplementary_Information_unformatted.docx \
"Figure S1%6%results/suppfigures/FigureS1.sample_UMAP.jpg" \
"Figure S2%6.5%results/suppfigures/FigureS2.SLITRK5_IGV.1.png" \
"Figure S2%6.5%results/suppfigures/FigureS2.SLITRK5_IGV.2.png" \
"Figure S2%6.5%results/suppfigures/FigureS2.SLITRK5_IGV.3.png" \
"Figure S2%6.5%results/suppfigures/FigureS2.SLITRK5_IGV.4.png" \
"Figure S2%6.5%results/suppfigures/FigureS2.SLITRK5_IGV.5.png" \
"Figure S2%6.5%results/suppfigures/FigureS2.SLITRK5_IGV.6.png" \
"Figure S3%6%results/suppfigures/FigureS3.SLITRK5_lollipop_plot.png" \
"Figure S4%6%results/suppfigures/FigureS4.clps_caco_rate_pca_PC1-PC8.jpg" \
"Figure S5%6%results/suppfigures/FigureS5.LOF_LOEUF_linear_fet.png" \
"Figure S6%6%results/suppfigures/FigureS6.LOF_PLI_LogisticLinearFet.png" \
"Figure S7%6%results/suppfigures/FigureS7.DNM_caco_rate.png" \
"Figure S8%6%results/suppfigures/FigureS8.DNM_rate.png" \
"Figure S9%6%results/suppfigures/FigureS9.DNM_rate.non_splice_frameshift.LOFintoltol.png" \
"Figure S10%6%results/suppfigures/FigureS10.DNM_misDLOF_rate_genesets.png" \
"Figure S11%6%results/suppfigures/FigureS11.DNM_misD_LOEUF_bins.png" \
"Figure S12%6%results/suppfigures/FigureS12.misD_MPC_bins_caco.png" \
"Figure S13%6%results/suppfigures/FigureS13.DNM_LOFloeuf1_misDmpcgt2_carriers_features_assoc.png" \
"Figure S14%6%results/suppfigures/FigureS14.LoF_LOEUFdeciles_ORs.malevsfemale.nontrio.png" \
"Figure S15%5%results/suppfigures/FigureS15.CHD8_LoF_IGV_Sanger.png" \
"Figure S16%6%results/suppfigures/FigureS16.HDAC4_LoF_IGV_Sanger.png"

exit
