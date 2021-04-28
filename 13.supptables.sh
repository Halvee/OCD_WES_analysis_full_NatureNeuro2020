#!/bin/bash

# needed for some python scripts - access to installed packages
export PYTHONPATH="$PWD/lib/python/:$PYTHONPATH"

# needed for some R scripts - access to installed packages
export R_LIBS=$PWD/lib/R/

# make sure output dir exists
mkdir -p results/supptables/

# make full sampleped from all 11 case/control collapsings
cat data/clps_meta/cluster_*.caco.sampleped \
> results/supptables/caco.sampleped

# Table S1 : make the full cohort table with all sample info included
bin/Rscript src/supptables/compile_cohort_table.R \
results/supptables/caco.sampleped \
data/clps_rate/OCD_clps_rate.allcases.sampleped \
results/input/OCDfams_2019.cohort.triosquartets.sampleped \
data/misc/cohort_stats/ocd_fam_iids.in_analysis.txt \
data/supptables/all_ocd_samples.dragendb.20190812.csv \
results/supptables/TableS1.full_cohort_manifest.csv

# Table S2 : case/control information per analysis
rsync \
results/misc/cohort_stats/cohort_full_stats.info.csv \
results/supptables/TableS2.casecontrol_analysis_info.csv

# Table S3 : phenotype counts in control cohorts
rsync \
results/misc/cohort_stats/cohort_full_stats.broad_phenotype_counts.csv \
results/supptables/TableS3.control_phenotype_counts.csv

# Table S4 : collapsing meta-analysis results (LoF)
bin/python src/supptables/xlsx_to_csv.py \
data/clps_meta/All_CMH_exact_summary_lclust_res.xlsx \
LOF \
results/supptables/collapsing_meta_lof.unformatted.csv

bin/Rscript src/supptables/format_clps_meta_csv.R \
data/clps_meta/ \
results/supptables/collapsing_meta_lof.unformatted.csv \
results/input/Homo_sapiens.GRCh37.87.ensg_genesymbol_entrezid.tsv \
results/supptables/TableS4.collapsing_meta_lof.csv 

# Table S5 : collapsing meta-analysis results (LoF/MisD)
bin/python src/supptables/xlsx_to_csv.py \
data/clps_meta/All_CMH_exact_summary_lclust_res.xlsx \
FUNCTIONAL \
results/supptables/collapsing_meta_lofmisd.unformatted.csv

bin/Rscript src/supptables/format_clps_meta_csv.R \
data/clps_meta/ \
results/supptables/collapsing_meta_lofmisd.unformatted.csv \
results/input/Homo_sapiens.GRCh37.87.ensg_genesymbol_entrezid.tsv \
results/supptables/TableS5.collapsing_meta_lofmisd.csv 

# Table S6 : Qualifying variation in SLITRK5 (with missense constraint annots)
rsync \
results/misc/slitrk5_missense_constraint/ocd_collapsing.slitrk5_qv.mis_ann.csv \
results/supptables/TableS6.SLITRK5_qualifying_variants.csv

# Table S7 : burden of missense annotation in constrained regions
# (first testing baseline burden of ppn2 damaging, then testing
#  burden of ppn2 damaging + missense constraint > 75th percentile 
#  of SLITRK5 variants collected in cases and on controls
rsync \
results/misc/slitrk5_missense_constraint/slitrk5_missense_constraint.csv \
results/supptables/TableS7.SLITRK5_missense_constraint_casectrl.csv

# Table S12 : Sanger validation results for select SNV/indel calls
# Table S13 : Summary statistics for sanger validation results
bin/Rscript src/supptables/get_sanger_results_stats.R \
data/supptables/OCD_DNM_2019_sanger_results_full_amended.csv \
results/supptables/TableS12.sanger_validation_results.csv \
results/supptables/TableS13.sanger_validation_sumstats.csv

# Table S8, S9 : de novo SNV/indel calls from trios. quartets
bin/Rscript src/supptables/dnm_callset_to_trios_quartets.R \
results/input/OCDfams_2019.cohort.trios.sampleped \
results/input/OCDfams_2019.cohort.quartets.sampleped \
data/input/OCDfams_2019.dnm.callset_final.sanger_adjust.csv \
results/supptables/TableS12.sanger_validation_results.csv \
results/dnm_casectrl/OCDfams_2019.cohort.trios.10x.ppn2hdiv.inkits.cdnm \
results/dnm_casectrl/OCDfams_2019.cohort.quartets.10x.ppn2hdiv.inkits.cdnm \
results/input/Homo_sapiens.GRCh37.87.ensg_genesymbol_entrezid.tsv \
results/supptables/TableS8.denovo_snvindel_calls.trios.csv \
results/supptables/TableS9.denovo_snvindel_calls.quartets.csv

# Table S10 : trio joint coverage summary by sample
bin/Rscript src/supptables/tsv_to_csv.R \
results/trio_coverage/OCDfams_2019.trios.cvg_stats.10x.sampleProfileSummary.tsv \
results/supptables/TableS10.10x_joint_cvg_per_trio.csv

# Table S11 trio joint coverage summary by CCDSr20 gene
bin/Rscript src/supptables/tsv_to_csv.R \
results/trio_coverage/OCDfams_2019.trios.cvg_stats.10x.genicProfileSummary.tsv \
results/supptables/TableS11.10x_joint_cvg_per_CCDSr20_gene.csv

# Table S14 : Full set of annotations from JHU case trios, 
#             JHU case quartets, Iossifov control trios, 
#             and additional Cappi case trios
bin/Rscript src/supptables/cdnms_to_csv.R \
results/dnm_casectrl/OCDfams_2019.cohort.trios.10x.ppn2hdiv.inkits.cdnm \
results/dnm_casectrl/OCDfams_2019.cohort.quartets.10x.ppn2hdiv.inkits.cdnm \
results/dnm_casectrl/Iossifov_etal_2014_healthysib_DNMs.ppn2hdiv.inkits.cdnm \
results/dnm_casectrl/Cappi_etal_2019_OCD.DNMs.ppn2hdiv.inkits.cdnm \
results/input/Homo_sapiens.GRCh37.87.ensg_genesymbol_entrezid.tsv \
results/supptables/TableS14.dnm_analysis_calls.csv

# Table S15 : DNM gene-based test results + extTADA results
bin/Rscript src/supptables/merge_dnm_test_results.R \
results/input/Homo_sapiens.GRCh37.87.ensg_genesymbol_entrezid.tsv \
results/dnm_genebased/OCD_JHU_Cappi_n771.lof.all.tsv \
results/dnm_genebased/OCD_JHU_Cappi_n771.dmg.all.tsv \
results/extTADA/OCD.771trios_476ca_1761co.extTADA.input.tsv \
results/extTADA/OCD.771trios_476ca_1761co.final.extTADA.20000perm.gene_res.tsv \
results/supptables/TableS15.dnmgenebased_extTADA_results.csv

# Table S16 : highlights of the more notable DNMs detected
rsync \
results/misc/notable_dnms/OCD.trios_quartets.notable_DNMs.csv \
results/supptables/TableS16.notable_DNMs.csv

# combine all supp tables into one xlsx file for distribution
bin/python src/supptables/csvs_to_xlsx.py \
--legend-txt data/supptables/tables.legend.txt \
--out-xlsx results/supptables/Supplemental_Tables_unformatted.xlsx \
"Table S1%results/supptables/TableS1.full_cohort_manifest.csv" \
"Table S2%results/supptables/TableS2.casecontrol_analysis_info.csv" \
"Table S3%results/supptables/TableS3.control_phenotype_counts.csv" \
"Table S4%results/supptables/TableS4.collapsing_meta_lof.csv" \
"Table S5%results/supptables/TableS5.collapsing_meta_lofmisd.csv" \
"Table S6%results/supptables/TableS6.SLITRK5_qualifying_variants.csv" \
"Table S7%results/supptables/TableS7.SLITRK5_missense_constraint_casectrl.csv" \
"Table S8%results/supptables/TableS8.denovo_snvindel_calls.trios.csv" \
"Table S9%results/supptables/TableS9.denovo_snvindel_calls.quartets.csv" \
"Table S10%results/supptables/TableS10.10x_joint_cvg_per_trio.csv" \
"Table S11%results/supptables/TableS11.10x_joint_cvg_per_CCDSr20_gene.csv" \
"Table S12%results/supptables/TableS12.sanger_validation_results.csv" \
"Table S13%results/supptables/TableS13.sanger_validation_sumstats.csv" \
"Table S14%results/supptables/TableS14.dnm_analysis_calls.csv" \
"Table S15%results/supptables/TableS15.dnmgenebased_extTADA_results.csv" \
"Table S16%results/supptables/TableS16.notable_DNMs.csv"

exit
