#!/bin/bash

# make sure output dir exists
mkdir -p results/mainfigures/

# figure 2 : collapsing meta-analysis results
# 2A) LOF dominant model
# 2B) LOF/misD dominant model
cd data/mainfigures/
rsync -L \
All_LOF_CMH_qq_exact_res_0_4_min_sample_5.png \
../../results/mainfigures/Figure2A.clps_meta.LOF.png
rsync -L \
All_FUNCTIONAL_CMH_qq_exact_res_0_4_min_sample_5.png \
../../results/mainfigures/Figure2B.clps_meta.LOFmisD.png
cd ../../

# figure 3 : case/control rate comparison of
# LOF variants in LOEUF bins (odds ratios)
rsync \
results/clps_rate/OCD_clps_rate.allcases.LOF_LOEUF.logistic.pdf \
results/mainfigures/Figure3.clps_rate.LOF_LOEUF.logistic.pdf

# figure 4 : case/control de novo mutation burden across
# the exome by annotation (top), for LOF variants within
# LOEUF bins (bottom)
rsync \
results/dnm_rate_caco/OCDfams_2019.cohort.trios.10x.caco.pdf \
results/mainfigures/Figure4.DNM_burden_vs_ctrls.logistic.pdf

# figure 5 : gene-based test results for DNM
bin/Rscript src/mainfigures/plot_extTADA.R  \
results/input/addjusted.CCDS.genes.index.r20.hg19.bed \
results/extTADA/OCD.771trios_476ca_1761co.final.extTADA.20000perm.gene_res.tsv \
results/mainfigures/Figure5.extTADA.manhattan.pdf

exit
