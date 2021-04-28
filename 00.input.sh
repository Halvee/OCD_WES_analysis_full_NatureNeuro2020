#!/bin/bash

# needed for some python scripts - access to installed packages
export PYTHONPATH="$PWD/lib/python/:$PYTHONPATH"

# needed for some R scripts - access to installed packages
export R_LIBS=$PWD/lib/R/

# make sure output dirs exist
mkdir -p results/input/
mkdir -p results/clps_rate/
mkdir -p results/mut_rate_tbls/
mkdir -p results/trio_coverage/
mkdir -p results/dnm_casectrl/
mkdir -p results/dnm_rate_caco/
mkdir -p results/dnm_genebased/
mkdir -p results/extTADA/

# create a masterset of gene symbols from CCDSr20
awk '{print $1}' data/input/addjusted.CCDS.genes.index.r20.hg19.txt \
| sort \
| uniq \
> results/input/CCDSr20.geneset

# get subset of CCDSr20 gene symbols that are on autosomes
awk '{if ($2!="X" && $2!="Y") {print $1}}' \
  data/input/addjusted.CCDS.genes.index.r20.hg19.txt \
| sort \
| uniq \
> results/input/CCDSr20_autosomal.geneset

# make a bed file with intervals from CCDSr20 masterfile
bin/python src/input/rgn_to_gene_bed.py \
data/input/addjusted.CCDS.genes.index.r20.hg19.txt \
| sort -k1,1 -k2,2n \
> results/input/addjusted.CCDS.genes.index.r20.hg19.bed

# download gnomAD LOEUF by gene
cd results/input/
wget -N https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
zcat gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
| awk -F"\t" '{OFS="\t"; print $1,$64,$36}' \
> gnomad.v2.1.1.loeuf_bins.tsv
zcat gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
| awk -F"\t" '{OFS="\t"; print $1,$21}' \
> gnomad.v2.1.1.pLI.tsv
cd ../../

# convert from xlsx to tsv (Iossifov et al. 2014)
bin/python src/input/format_Iossifov_2014_DNMs.py \
data/input/NIHMS659757-supplement-Supplementary_Table_2.xlsx \
results/input/Iossifov_etal_2014_healthysib_DNMs.tsv

# convert from xlsx to tsv (Cappi et al. 2019)
bin/python src/input/format_Cappi_2019_DNMs.py \
data/input/TableS2_DN\ variants.xlsx \
results/input/Cappi_etal_2019_OCD.DNMs.tsv

# convert from xlsx to tsv (Iossifov et al. 2014, shared DNMs btwn sibs specifically)
bin/python src/input/format_Iossifov_2014_sharedDNMs.py \
data/input/NIHMS659757-supplement-Supplementary_Table_2.xlsx \
results/input/Iossifov_etal_2014_sibling_sharedDNMs.tsv

# get varlists for MPC 0-1, 1-2, >2
zcat data/input/fordist_constraint_official_only_mpc_values_v2.tsv.gz \
| awk '{ if ((NR!=1) && ($5<=1) && ($5!="NA")){print $1"-"$2"-"$3"-"$4}}' \
| gzip -c - \
> results/input/MPC0to1.varlist.gz
zcat data/input/fordist_constraint_official_only_mpc_values_v2.tsv.gz \
| awk '{ if ((NR!=1) && ($5>1) && ($5<=2) && ($5!="NA")){print $1"-"$2"-"$3"-"$4}}' \
| gzip -c - \
> results/input/MPC1to2.varlist.gz
zcat data/input/fordist_constraint_official_only_mpc_values_v2.tsv.gz \
| awk '{ if ((NR!=1) && ($5>2) && ($5!="NA")){print $1"-"$2"-"$3"-"$4}}' \
| gzip -c - \
> results/input/MPCgt2.varlist.gz

# convert raw table contents to geneset file (n=124 genes that
# are exomewide significant, ie. p<5e-7)
cat data/input/PMID_30559488_Table2.raw.txt \
| cut -f2 \
| tr "," "\n" \
| sed 's/\ *//g' \
| sort \
| uniq \
| grep "*$" \
| sed 's/\*$//g; s/a$//g' \
> results/input/PMID_30559488_Table2.geneset

# get set of Neurodevelopmental disorder genes that fit one of the criteria
# listed here:
# 1. a gene with q < 0.1 from Satterstrom et al. 2020 ASD DNM study 
#    (n=102 genes)
# 2. a gene listed as exomewide significant in Coe et al. 2019 for assoc
#    with neurodevelopmental disease
#    (n=124 genes)
# Total number of genes : 187
bin/python \
src/input/form_NDD_geneset_file.py \
data/input/1-s2.0-S0092867419313984-mmc2.xlsx \
results/input/PMID_30559488_Table2.geneset \
results/input/gnomad.v2.1.1.loeuf_bins.tsv \
results/input/NDD.n187.csv \
results/input/NDD.n187.geneset \
> results/input/NDD.n187.log

# read gene-based results from Wang et al. 2018 study of ~ 800 TS trios, make
# two seperate geneset files:
# 1. genes hit with 2 or more damaging de novos (misD, LOF)
# 2. genes hit with 1 damaging de novo (misD, LOF)
bin/python \
src/input/get_Wang_2018_TS_genes.py \
data/input/1-s2.0-S221112471831386X-mmc4.xlsx \
results/input/TS.ge2_dmg_DNMs_q_lt30.geneset \
results/input/TS.ge1_dmg_DNM_q_ge30.geneset \
> results/input/TS.log

# assemble input trio and quartet analysis cohort
bin/Rscript src/input/dnm_cohort_prepare.R \
data/input/OCDfams_2019.triosquartets.analysisready.sampleped \
data/input/OCDfams_2019.dnm.callset_final.sanger_adjust.csv \
results/input/OCDfams_2019.cohort \
> results/input/OCDfams_2019.cohort.log

# make table of gene symbols / entrez ids / ensembl gene ids
#bin/Rscript src/input/gene_symbol_entrez_ensg_table.R \
# make table of ensembl gene id / gene symbol combinations from Ensembl GRCh 37.87
zcat data/input/Homo_sapiens.GRCh37.87.gtf.gz \
| grep ensembl \
| grep protein_coding \
| awk '{if ($3=="gene") {print $10,$14}}' \
| tr -d '"' \
| tr -d ";" \
| sort \
| uniq \
> results/input/Homo_sapiens.GRCh37.87.ensg_genesymbol.tsv

# use table as input to get additional entrez IDs, return a table with only unique
# ensembl gene id / gene symbol / entrez id values
bin/Rscript src/input/gene_symbol_entrez_ensg_table.R \
results/input/Homo_sapiens.GRCh37.87.ensg_genesymbol.tsv \
results/input/CCDSr20.geneset \
results/input/Homo_sapiens.GRCh37.87.ensg_genesymbol_entrezid.tsv

exit
